// ShearPrjGPU.cu -- GPU module for Costanzi-2026 shear projection integrals
//
// Computes sigma_prj, dsigma_prj, and shear_prj (gamma_t^prj) on GPU.
// Uses PAGANI GPU integrator like bSelMargGPU for proper device execution.
//
// Outputs to datablock:
//   sigma_prj/{vals, rnd, cl}
//   dsigma_prj/{vals, rnd, cl}
//   shear_prj/{vals, rnd, cl}

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"

#include "utils/make_grid_points.hh"
#include "utils/datablock_reader.hh"
#include "utils/gpu_integrator.cuh"
#include "utils/mpi_support.hh"

#include "models/sigma_prj_gpu_t.cuh"

#include <iostream>
#include <optional>
#include <vector>

using cosmosis::DataBlock;

namespace {

  struct ShearPrjGPUConfig {
    quad::Volume<double, 3> volume;
    y3_cluster::grid_t<4> grid_points;
    double eps_rel;
    double eps_abs;
    double max_eval;
    std::string algorithm;
    int device_id;
  };

  ShearPrjGPUConfig
  read_config(DataBlock& cfg)
  {
    char const* label = "ShearPrjGPU";
    ShearPrjGPUConfig c;

    // Integration bounds: (theta, z, lnM)
    c.volume.lows[0]  = cfg.view<double>(label, "theta_low");
    c.volume.highs[0] = cfg.view<double>(label, "theta_high");
    c.volume.lows[1]  = cfg.view<double>(label, "zt_low");
    c.volume.highs[1] = cfg.view<double>(label, "zt_high");
    c.volume.lows[2]  = cfg.view<double>(label, "lnm_low");
    c.volume.highs[2] = cfg.view<double>(label, "lnm_high");

    // Grid points: (lambda_bin, zo_low, zo_high, radii) zipped
    auto lambda_bin = get_vector_double(cfg, label, "lambda_bin");
    auto zo_low     = get_vector_double(cfg, label, "zo_low");
    auto zo_high    = get_vector_double(cfg, label, "zo_high");
    auto radii      = get_vector_double(cfg, label, "radii");

    c.grid_points.set_names({"lambda_bin", "zo_low", "zo_high", "radii"});
    for (size_t i = 0; i < lambda_bin.size(); ++i) {
      c.grid_points.push_back({lambda_bin[i], zo_low[i], zo_high[i], radii[i]});
    }

    // Integration settings
    c.eps_rel   = cfg.view<double>(label, "eps_rel");
    c.eps_abs   = cfg.view<double>(label, "eps_abs");
    c.max_eval  = cfg.view<double>(label, "max_eval");
    c.algorithm = cfg.view<std::string>(label, "algorithm");

    // CUDA device
    auto info = y3_cluster::get_mpi_info();
    c.device_id = info.local_rank;

    return c;
  }

}  // anonymous namespace


class ShearPrjGPU {
public:
  explicit ShearPrjGPU(DataBlock& cfg)
    : cfg_(read_config(cfg))
  {
    cudaSetDevice(cfg_.device_id);
  }

  void execute(DataBlock& sample)
  {
    cudaSetDevice(cfg_.device_id);

    size_t const ngrid = cfg_.grid_points.size();

    // Storage for results
    std::vector<double> sigma_rnd(ngrid), sigma_cl(ngrid), sigma_vals(ngrid);
    std::vector<double> dsigma_rnd(ngrid), dsigma_cl(ngrid), dsigma_vals(ngrid);
    std::vector<double> shear_rnd(ngrid), shear_cl(ngrid), shear_vals(ngrid);

    // Create GPU integrator
    y3_cuda::MultiDimensionalIntegrator<3> integrator(cfg_.algorithm);
    integrator.set_maxeval(cfg_.max_eval);

    // Create integrands for each component
    y3_cuda::SigmaPrj_gpu<y3_cuda::SigmaRndWeight>  sig_rnd_integrand(sample);
    y3_cuda::SigmaPrj_gpu<y3_cuda::SigmaClWeight>   sig_cl_integrand(sample);
    y3_cuda::SigmaPrj_gpu<y3_cuda::DSigmaRndWeight> dsig_rnd_integrand(sample);
    y3_cuda::SigmaPrj_gpu<y3_cuda::DSigmaClWeight>  dsig_cl_integrand(sample);

    // Initialize models from sample
    sig_rnd_integrand.set_sample(sample);
    sig_cl_integrand.set_sample(sample);
    dsig_rnd_integrand.set_sample(sample);
    dsig_cl_integrand.set_sample(sample);

    // Integrate for each grid point
    for (size_t ig = 0; ig < ngrid; ++ig) {
      auto const& gp = cfg_.grid_points.points[ig];

      // Sigma random
      sig_rnd_integrand.set_grid_point(gp);
      auto res_sr = integrator.integrate(sig_rnd_integrand, cfg_.eps_abs,
                                          cfg_.eps_rel, cfg_.volume);
      sigma_rnd[ig] = res_sr.estimate;

      // Sigma clustering
      sig_cl_integrand.set_grid_point(gp);
      auto res_sc = integrator.integrate(sig_cl_integrand, cfg_.eps_abs,
                                          cfg_.eps_rel, cfg_.volume);
      sigma_cl[ig] = res_sc.estimate;

      // DSigma random
      dsig_rnd_integrand.set_grid_point(gp);
      auto res_dr = integrator.integrate(dsig_rnd_integrand, cfg_.eps_abs,
                                          cfg_.eps_rel, cfg_.volume);
      dsigma_rnd[ig] = res_dr.estimate;

      // DSigma clustering
      dsig_cl_integrand.set_grid_point(gp);
      auto res_dc = integrator.integrate(dsig_cl_integrand, cfg_.eps_abs,
                                          cfg_.eps_rel, cfg_.volume);
      dsigma_cl[ig] = res_dc.estimate;

      // Totals
      sigma_vals[ig]  = sigma_rnd[ig] + sigma_cl[ig];
      dsigma_vals[ig] = dsigma_rnd[ig] + dsigma_cl[ig];

      // Shear = DSigma * Sigma_crit_inv
      double const sci = sig_rnd_integrand.get_sci();
      shear_rnd[ig]  = dsigma_rnd[ig] * sci;
      shear_cl[ig]   = dsigma_cl[ig] * sci;
      shear_vals[ig] = dsigma_vals[ig] * sci;
    }

    // Write results to datablock
    using cosmosis::ndarray;

    // sigma_prj
    ndarray<double> sig_v(sigma_vals, {ngrid});
    ndarray<double> sig_r(sigma_rnd, {ngrid});
    ndarray<double> sig_c(sigma_cl, {ngrid});
    sample.put_val("sigma_prj", "vals", sig_v);
    sample.put_val("sigma_prj", "rnd", sig_r);
    sample.put_val("sigma_prj", "cl", sig_c);

    // dsigma_prj
    ndarray<double> dsig_v(dsigma_vals, {ngrid});
    ndarray<double> dsig_r(dsigma_rnd, {ngrid});
    ndarray<double> dsig_c(dsigma_cl, {ngrid});
    sample.put_val("dsigma_prj", "vals", dsig_v);
    sample.put_val("dsigma_prj", "rnd", dsig_r);
    sample.put_val("dsigma_prj", "cl", dsig_c);

    // shear_prj (gamma_t^prj)
    ndarray<double> sh_v(shear_vals, {ngrid});
    ndarray<double> sh_r(shear_rnd, {ngrid});
    ndarray<double> sh_c(shear_cl, {ngrid});
    sample.put_val("shear_prj", "vals", sh_v);
    sample.put_val("shear_prj", "rnd", sh_r);
    sample.put_val("shear_prj", "cl", sh_c);
  }

private:
  ShearPrjGPUConfig cfg_;
};


// CosmoSIS interface
extern "C" {

void* setup(DataBlock* cfg)
{
  return new ShearPrjGPU(*cfg);
}

DATABLOCK_STATUS execute(DataBlock* sample, void* module)
{
  auto mod = static_cast<ShearPrjGPU*>(module);
  mod->execute(*sample);
  return DBS_SUCCESS;
}

int cleanup(void* module)
{
  delete static_cast<ShearPrjGPU*>(module);
  return 0;
}

}  // extern "C"
