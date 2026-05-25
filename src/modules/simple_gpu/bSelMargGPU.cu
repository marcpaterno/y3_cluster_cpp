// bSelMargGPU.cu -- Unified GPU module for Costanzi-2026 (P1, I1, J) operators
//
// Adapts p_operator_cuhre_t.hh for GPU execution. Computes all three scalars
// in a single module to share model setup overhead.
//
// Outputs to datablock:
//   bSelMargGPU/P1_vals, P1_errors
//   bSelMargGPU/I1_vals, I1_errors
//   bSelMargGPU/J_vals,  J_errors
//   bSelMargGPU/zo_low, zo_high, lambda_bin  (grid points)

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"

#include "utils/make_grid_points.hh"
#include "utils/datablock_reader.hh"
#include "utils/gpu_integrator.cuh"
#include "utils/mpi_support.hh"

#include "models/p_operator_gpu_t.cuh"

#include <iostream>
#include <optional>
#include <vector>

using cosmosis::DataBlock;

namespace {

  // Module configuration
  struct BSelMargGPUConfig {
    quad::Volume<double, 3> volume;
    y3_cluster::grid_t<3> grid_points;
    double eps_rel;
    double eps_abs;
    double max_eval;
    std::string algorithm;
    int device_id;
  };

  BSelMargGPUConfig
  read_config(DataBlock& cfg)
  {
    char const* label = "bSelMargGPU";
    BSelMargGPUConfig c;

    // Integration bounds
    c.volume.lows[0]  = cfg.view<double>(label, "lt_low");
    c.volume.highs[0] = cfg.view<double>(label, "lt_high");
    c.volume.lows[1]  = cfg.view<double>(label, "zt_low");
    c.volume.highs[1] = cfg.view<double>(label, "zt_high");
    c.volume.lows[2]  = cfg.view<double>(label, "lnm_low");
    c.volume.highs[2] = cfg.view<double>(label, "lnm_high");

    // Grid points
    auto zo_low_vec  = get_vector_double(cfg, label, "zo_low");
    auto zo_high_vec = get_vector_double(cfg, label, "zo_high");
    auto bin_vec     = get_vector_double(cfg, label, "lambda_bin");

    c.grid_points.set_names({"zo_low", "zo_high", "lambda_bin"});
    for (size_t i = 0; i < zo_low_vec.size(); ++i) {
      c.grid_points.push_back({zo_low_vec[i], zo_high_vec[i], bin_vec[i]});
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


class BSelMargGPU {
public:
  explicit BSelMargGPU(DataBlock& cfg)
    : cfg_(read_config(cfg))
  {
    cudaSetDevice(cfg_.device_id);
  }

  void execute(DataBlock& sample)
  {
    cudaSetDevice(cfg_.device_id);

    size_t const ngrid = cfg_.grid_points.size();

    // Storage for results
    std::vector<double> P1_vals(ngrid), P1_errs(ngrid);
    std::vector<double> I1_vals(ngrid), I1_errs(ngrid);
    std::vector<double> J_vals(ngrid),  J_errs(ngrid);

    // Create integrator
    y3_cuda::MultiDimensionalIntegrator<3> integrator(cfg_.algorithm);
    integrator.set_maxeval(cfg_.max_eval);

    // Create integrands for each operator type
    // They share the same underlying models but use different weight functions
    y3_cuda::P_operator_gpu<y3_cuda::BSelMargP1GPUWeight> p1_integrand(sample);
    y3_cuda::P_operator_gpu<y3_cuda::BSelMargI1GPUWeight> i1_integrand(sample);
    y3_cuda::P_operator_gpu<y3_cuda::BSelMargJGPUWeight>  j_integrand(sample);

    // Initialize models from sample (done once, shared via sample)
    p1_integrand.set_sample(sample);
    i1_integrand.set_sample(sample);
    j_integrand.set_sample(sample);

    // Integrate for each grid point
    for (size_t ig = 0; ig < ngrid; ++ig) {
      auto const& gp = cfg_.grid_points.points[ig];

      // P1
      p1_integrand.set_grid_point(gp);
      auto res_p1 = integrator.integrate(p1_integrand, cfg_.eps_abs, cfg_.eps_rel, cfg_.volume);
      P1_vals[ig] = res_p1.estimate;
      P1_errs[ig] = res_p1.errorest;

      // I1
      i1_integrand.set_grid_point(gp);
      auto res_i1 = integrator.integrate(i1_integrand, cfg_.eps_abs, cfg_.eps_rel, cfg_.volume);
      I1_vals[ig] = res_i1.estimate;
      I1_errs[ig] = res_i1.errorest;

      // J = I2 - I1 (computed directly)
      j_integrand.set_grid_point(gp);
      auto res_j = integrator.integrate(j_integrand, cfg_.eps_abs, cfg_.eps_rel, cfg_.volume);
      J_vals[ig] = res_j.estimate;
      J_errs[ig] = res_j.errorest;
    }

    // Write results to datablock
    using cosmosis::ndarray;
    char const* section = "bSelMargGPU";

    ndarray<double> P1_arr(P1_vals, {ngrid});
    ndarray<double> P1_err_arr(P1_errs, {ngrid});
    sample.put_val(section, "P1_vals", P1_arr);
    sample.put_val(section, "P1_errors", P1_err_arr);

    ndarray<double> I1_arr(I1_vals, {ngrid});
    ndarray<double> I1_err_arr(I1_errs, {ngrid});
    sample.put_val(section, "I1_vals", I1_arr);
    sample.put_val(section, "I1_errors", I1_err_arr);

    ndarray<double> J_arr(J_vals, {ngrid});
    ndarray<double> J_err_arr(J_errs, {ngrid});
    sample.put_val(section, "J_vals", J_arr);
    sample.put_val(section, "J_errors", J_err_arr);

    // Write grid points for downstream use
    std::vector<double> zo_low_out(ngrid), zo_high_out(ngrid), bin_out(ngrid);
    for (size_t i = 0; i < ngrid; ++i) {
      zo_low_out[i]  = cfg_.grid_points.points[i][0];
      zo_high_out[i] = cfg_.grid_points.points[i][1];
      bin_out[i]     = cfg_.grid_points.points[i][2];
    }
    ndarray<double> zo_low_arr(zo_low_out, {ngrid});
    ndarray<double> zo_high_arr(zo_high_out, {ngrid});
    ndarray<double> bin_arr(bin_out, {ngrid});
    sample.put_val(section, "zo_low", zo_low_arr);
    sample.put_val(section, "zo_high", zo_high_arr);
    sample.put_val(section, "lambda_bin", bin_arr);
  }

private:
  BSelMargGPUConfig cfg_;
};


// CosmoSIS interface
extern "C" {

void* setup(DataBlock* cfg)
{
  return new BSelMargGPU(*cfg);
}

DATABLOCK_STATUS execute(DataBlock* sample, void* module)
{
  auto mod = static_cast<BSelMargGPU*>(module);
  mod->execute(*sample);
  return DBS_SUCCESS;
}

int cleanup(void* module)
{
  delete static_cast<BSelMargGPU*>(module);
  return 0;
}

}  // extern "C"
