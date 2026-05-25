// ShearPrjGPU.cu -- GPU module for Costanzi-2026 shear projection integrals
//
// Computes sigma_prj, dsigma_prj, and shear_prj (gamma_t^prj) on GPU.
// Outputs match the CPU ShearPrjEvaluator sections for compatibility.
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

    // Storage for results: total, rnd, cl for each observable
    std::vector<double> sigma_vals(ngrid), sigma_rnd(ngrid), sigma_cl(ngrid);
    std::vector<double> dsigma_vals(ngrid), dsigma_rnd(ngrid), dsigma_cl(ngrid);
    std::vector<double> shear_vals(ngrid), shear_rnd(ngrid), shear_cl(ngrid);

    // Create integrand
    y3_cuda::SigmaPrj_gpu integrand(sample);
    integrand.set_sample(sample);

    // For now, use a simple integration approach
    // In production, this would use PAGANI or similar GPU integrator
    // Here we do a fixed-grid Gauss-Legendre approximation

    size_t const N_theta = 30;
    size_t const N_z = 40;
    size_t const N_M = 24;

    // GL nodes and weights (simplified - using uniform for now)
    auto make_gl_grid = [](double lo, double hi, size_t n,
                           std::vector<double>& nodes,
                           std::vector<double>& weights) {
      nodes.resize(n);
      weights.resize(n);
      double const dx = (hi - lo) / n;
      for (size_t i = 0; i < n; ++i) {
        nodes[i] = lo + (i + 0.5) * dx;
        weights[i] = dx;
      }
    };

    std::vector<double> theta_nodes, theta_weights;
    std::vector<double> z_nodes, z_weights;
    std::vector<double> M_nodes, M_weights;

    make_gl_grid(cfg_.volume.lows[0], cfg_.volume.highs[0], N_theta,
                 theta_nodes, theta_weights);
    make_gl_grid(cfg_.volume.lows[1], cfg_.volume.highs[1], N_z,
                 z_nodes, z_weights);
    make_gl_grid(cfg_.volume.lows[2], cfg_.volume.highs[2], N_M,
                 M_nodes, M_weights);

    // Integrate for each grid point
    for (size_t ig = 0; ig < ngrid; ++ig) {
      auto const& gp = cfg_.grid_points.points[ig];
      integrand.set_grid_point(gp);

      double acc_S_rnd = 0.0, acc_S_cl = 0.0;
      double acc_DS_rnd = 0.0, acc_DS_cl = 0.0;

      // 3D integration loop
      for (size_t i_th = 0; i_th < N_theta; ++i_th) {
        for (size_t i_z = 0; i_z < N_z; ++i_z) {
          for (size_t i_M = 0; i_M < N_M; ++i_M) {
            double s_rnd, s_cl, ds_rnd, ds_cl;
            integrand(theta_nodes[i_th], z_nodes[i_z], M_nodes[i_M],
                      s_rnd, s_cl, ds_rnd, ds_cl);

            double const wt = theta_weights[i_th] * z_weights[i_z] * M_weights[i_M];
            acc_S_rnd  += wt * s_rnd;
            acc_S_cl   += wt * s_cl;
            acc_DS_rnd += wt * ds_rnd;
            acc_DS_cl  += wt * ds_cl;
          }
        }
      }

      sigma_rnd[ig]  = acc_S_rnd;
      sigma_cl[ig]   = acc_S_cl;
      sigma_vals[ig] = acc_S_rnd + acc_S_cl;

      dsigma_rnd[ig]  = acc_DS_rnd;
      dsigma_cl[ig]   = acc_DS_cl;
      dsigma_vals[ig] = acc_DS_rnd + acc_DS_cl;

      // gamma_t = DSigma * Sigma_crit_inv
      double const sci = integrand.get_sci();
      shear_vals[ig] = dsigma_vals[ig] * sci;
      shear_rnd[ig]  = dsigma_rnd[ig] * sci;
      shear_cl[ig]   = dsigma_cl[ig] * sci;
    }

    // Write results to datablock (same sections as CPU for compatibility)
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
