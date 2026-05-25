// Total surface density (1h + b*2h) GPU module WITHOUT selection function.
//
// This module computes the bin-averaged total surface density:
//   <Sigma_tot>(R) = N^{-1} * int dlnM int dz  Omega(z) * (dV/dOmega/dz)(z)
//                    * n(M,z) * P(lo|M,z) * K(zo|z) * Sigma_tot(R, M, z)
//
// where Sigma_tot = Sigma_1h(R,M) + b(M,z) * Sigma_2h(R,z)
//
// Divide output by numberCountsSimple to get the bin-averaged value.
#include "utils/cuda_module_macros.cuh"
#include "utils/datablock_reader.hh"
#include "utils/make_cuda_integration_volumes.cuh"
#include "utils/make_grid_points.hh"

#include "cubacpp/integration_result.hh"
#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/datablock_status.h"
#include "common/cuda/Volume.cuh"

#include "models/omega_z_des.cuh"
#include "models/dv_do_dz_t.cuh"
#include "models/hmf_t.cuh"
#include "models/mor_des_log_t.cuh"
#include "models/int_zo_zt_des_t.cuh"
#include "models/sig_max.cuh"

#include <iostream>
#include <optional>
#include <vector>

using cosmosis::DataBlock;
using cubacpp::integration_result;

// sigmaTotSimple: Total surface density without selection function.
// Integrates over (lo, zt, lnM) with MOR and Sigma_tot(R, M, z).
class sigmaTotSimple {
public:
  using grid_t = y3_cluster::grid_t<3>;
  using grid_point_t = grid_t::value_type;

private:
  using volume_t = quad::Volume<double, 3>;

  // Models from datablock
  std::optional<y3_cuda::OMEGA_Z_DES> omega_z_;
  std::optional<y3_cuda::DV_DO_DZ_t> dv_do_dz_;
  std::optional<y3_cuda::HMF_t> hmf_;
  std::optional<y3_cuda::MOR_DES_LOG_t> mor_;
  std::optional<y3_cuda::INT_ZO_ZT_DES_t> int_zo_zt_;
  std::optional<y3_cuda::SIG_MAX> sigma_;

  // Current bin: (zo_low, zo_high, radius)
  double zo_low_;
  double zo_high_;
  double radius_;

  bool do_cartesian_product_of_bins_ = false;

public:
  explicit sigmaTotSimple(cosmosis::DataBlock& cfg)
  {
    auto rc = cfg.get_val(module_label(),
                          "do_cartesian_product_of_bins",
                          false,
                          do_cartesian_product_of_bins_);
    if (rc != DBS_SUCCESS) {
      throw std::runtime_error(
        "sigmaTotSimple: failed to read do_cartesian_product_of_bins\n");
    }
  }

  void set_sample(cosmosis::DataBlock& sample)
  {
    omega_z_.emplace(sample);
    dv_do_dz_.emplace(sample);
    hmf_.emplace(sample);
    mor_.emplace(sample);
    // INT_ZO_ZT_DES_t has default constructor (no DataBlock needed)
    int_zo_zt_.emplace();
    sigma_.emplace(sample);
  }

  void set_grid_point(grid_point_t const& grid_point)
  {
    zo_low_ = grid_point[0];
    zo_high_ = grid_point[1];
    radius_ = grid_point[2];
  }

  size_t get_device_mem_footprint()
  {
    size_t dev_size = 0;
    if (mor_) dev_size += mor_->get_device_mem_footprint();
    if (dv_do_dz_) dev_size += dv_do_dz_->get_device_mem_footprint();
    if (hmf_) dev_size += hmf_->get_device_mem_footprint();
    if (sigma_) dev_size += sigma_->get_device_mem_footprint();
    return dev_size;
  }

  __host__ __device__ double
  operator()(double lo, double zt, double lnM) const
  {
    // MOR: P(lambda_ob | M, z)
    double const mor_v = (*mor_)(lo, lnM, zt);

    // Common term: Omega(z) * dV/dOmega/dz * dn/dlnM
    double const common = (*omega_z_)(zt) * (*dv_do_dz_)(zt) * (*hmf_)(lnM, zt);

    // Photo-z kernel
    double const photoz = (*int_zo_zt_)(zo_low_, zo_high_, zt);

    // Total surface density: Sigma_1h + b*Sigma_2h
    double const sigma_v = (*sigma_)(radius_, lnM, zt);

    return common * mor_v * photoz * sigma_v;
  }

  static char const* module_label() { return "sigmaTotSimple"; }

  static std::vector<volume_t>
  make_integration_volumes(cosmosis::DataBlock& cfg)
  {
    return y3_cuda::make_integration_volumes_wall_of_numbers(
      cfg, sigmaTotSimple::module_label(), "lo", "zt", "lnm");
  }

  static grid_t make_grid_points(cosmosis::DataBlock& cfg)
  {
    return y3_cluster::make_grid_points_wall_of_numbers(
      cfg, sigmaTotSimple::module_label(), "zo_low", "zo_high", "radii");
  }
};

DEFINE_COSMOSIS_CUDA_INTEGRATION_MODULE(sigmaTotSimple)
