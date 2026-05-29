// Number counts GPU module with FULL 4D integration.
//
// This module computes cluster number counts by integrating over
// (lambda_ob, lambda_tr, z, lnM) directly WITHOUT using the pre-tabulated
// selection function S_ij.
//
// Integrand:
//   N = int dlob int dltr int dz int dlnM
//       Omega(z) * (dV/dOmega/dz)(z) * n(M,z)
//       * P_HOD(ltr|M,z) * P_EMG(lob|ltr,z) * K(zo|z)
//
// where:
//   - P_HOD(ltr|M,z) is the MOR (mass-observable relation): mor_des_log_t
//   - P_EMG(lob|ltr,z) is the EMG observation scatter: emg_des_t
//   - K(zo|z) is the photo-z bin kernel: int_zo_zt_des_t
//
// This 4D approach is used to validate the GPU pipeline independently
// from the CPU sel_function approach.
//
// Grid: (lob_low, lob_high, zo_low, zo_high) - per richness and redshift bin

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
#include "models/emg_des_t.cuh"

#include <iostream>
#include <optional>
#include <vector>

using cosmosis::DataBlock;
using cubacpp::integration_result;

// numberCountsFull: Full 4D integration over (lob, ltr, zt, lnM)
// without using pre-tabulated selection function.
class numberCountsFull {
public:
  // Grid point: (lob_low, lob_high, zo_low, zo_high)
  // This defines the richness bin [lob_low, lob_high] and redshift bin [zo_low, zo_high]
  using grid_t = y3_cluster::grid_t<4>;
  using grid_point_t = grid_t::value_type;

private:
  // 4D integration volume: (lob, ltr, zt, lnM)
  using volume_t = quad::Volume<double, 4>;

  // Models from datablock
  std::optional<y3_cuda::OMEGA_Z_DES> omega_z_;
  std::optional<y3_cuda::DV_DO_DZ_t> dv_do_dz_;
  std::optional<y3_cuda::HMF_t> hmf_;
  std::optional<y3_cuda::MOR_DES_LOG_t> mor_;        // P_HOD(ltr | M, z)
  std::optional<y3_cuda::EMG_DES_t> emg_;            // P_EMG(lob | ltr, z)
  std::optional<y3_cuda::INT_ZO_ZT_DES_t> int_zo_zt_;  // Photo-z kernel K(zo|z)

  // Current bin edges from grid point
  double lob_low_;
  double lob_high_;
  double zo_low_;
  double zo_high_;

public:
  explicit numberCountsFull(cosmosis::DataBlock& /*cfg*/)
  {
    // No special configuration needed at setup time
  }

  void set_sample(cosmosis::DataBlock& sample)
  {
    omega_z_.emplace(sample);
    dv_do_dz_.emplace(sample);
    hmf_.emplace(sample);
    mor_.emplace(sample);
    emg_.emplace(sample);
    // INT_ZO_ZT_DES_t has default constructor
    int_zo_zt_.emplace();
  }

  void set_grid_point(grid_point_t const& grid_point)
  {
    lob_low_ = grid_point[0];
    lob_high_ = grid_point[1];
    zo_low_ = grid_point[2];
    zo_high_ = grid_point[3];
  }

  size_t get_device_mem_footprint()
  {
    size_t dev_size = 0;
    if (omega_z_) dev_size += 0;  // OMEGA_Z_DES has no device memory
    if (dv_do_dz_) dev_size += dv_do_dz_->get_device_mem_footprint();
    if (hmf_) dev_size += hmf_->get_device_mem_footprint();
    if (mor_) dev_size += mor_->get_device_mem_footprint();
    if (emg_) dev_size += emg_->get_device_mem_footprint();
    return dev_size;
  }

  // The integrand: called at each (lob, ltr, zt, lnM) point
  // Integration variables are in the order expected by PAGANI
  __host__ __device__ double
  operator()(double lob, double ltr, double zt, double lnM) const
  {
    // Skip if lob is outside the current richness bin
    // (This allows using a single integration volume for multiple bins
    // by zeroing out contributions outside the bin)
    // Actually, the integration volume should already be set per-bin,
    // but we add this check for safety.
    if (lob < lob_low_ || lob > lob_high_) {
      return 0.0;
    }

    // Photo-z kernel: integral of P(z_ob | z_true) over [zo_low, zo_high]
    double const photoz = (*int_zo_zt_)(zo_low_, zo_high_, zt);
    if (photoz <= 0.0) {
      return 0.0;
    }

    // MOR: P_HOD(lambda_tr | M, z)
    double const p_hod = (*mor_)(ltr, lnM, zt);
    if (p_hod <= 0.0) {
      return 0.0;
    }

    // EMG: P(lambda_ob | lambda_tr, z)
    double const p_emg = (*emg_)(lob, ltr, zt);
    if (p_emg <= 0.0) {
      return 0.0;
    }

    // Common cosmology factors: Omega(z) * dV/dOmega/dz * dn/dlnM
    double const common = (*omega_z_)(zt) * (*dv_do_dz_)(zt) * (*hmf_)(lnM, zt);

    return common * p_hod * p_emg * photoz;
  }

  static char const* module_label() { return "numberCountsFull"; }

  // Create 4D integration volumes from ini config
  // The volume is over (lob, ltr, zt, lnM)
  static std::vector<volume_t>
  make_integration_volumes(cosmosis::DataBlock& cfg)
  {
    char const* label = module_label();

    // Read integration bounds
    double lob_low = cfg.view<double>(label, "lob_int_low");
    double lob_high = cfg.view<double>(label, "lob_int_high");
    double ltr_low = cfg.view<double>(label, "ltr_low");
    double ltr_high = cfg.view<double>(label, "ltr_high");
    double zt_low = cfg.view<double>(label, "zt_low");
    double zt_high = cfg.view<double>(label, "zt_high");
    double lnm_low = cfg.view<double>(label, "lnm_low");
    double lnm_high = cfg.view<double>(label, "lnm_high");

    // Single 4D volume for all bins
    // The actual bin limits are applied via grid points
    volume_t vol;
    vol.lows[0] = lob_low;
    vol.highs[0] = lob_high;
    vol.lows[1] = ltr_low;
    vol.highs[1] = ltr_high;
    vol.lows[2] = zt_low;
    vol.highs[2] = zt_high;
    vol.lows[3] = lnm_low;
    vol.highs[3] = lnm_high;

    return {vol};
  }

  // Create grid points from ini config
  // Each grid point is (lob_low, lob_high, zo_low, zo_high) for a bin
  static grid_t make_grid_points(cosmosis::DataBlock& cfg)
  {
    char const* label = module_label();

    // Read bin edges as wall-of-numbers
    auto lob_low_vec = get_vector_double(cfg, label, "lob_bin_low");
    auto lob_high_vec = get_vector_double(cfg, label, "lob_bin_high");
    auto zo_low_vec = get_vector_double(cfg, label, "zo_low");
    auto zo_high_vec = get_vector_double(cfg, label, "zo_high");

    // Check sizes match
    size_t n = lob_low_vec.size();
    if (lob_high_vec.size() != n || zo_low_vec.size() != n || zo_high_vec.size() != n) {
      throw std::runtime_error(
        "numberCountsFull: bin arrays must have same length");
    }

    grid_t result;
    result.set_names({"lob_bin_low", "lob_bin_high", "zo_low", "zo_high"});
    for (size_t i = 0; i < n; ++i) {
      result.push_back({lob_low_vec[i], lob_high_vec[i],
                        zo_low_vec[i], zo_high_vec[i]});
    }
    return result;
  }
};

DEFINE_COSMOSIS_CUDA_INTEGRATION_MODULE(numberCountsFull)
