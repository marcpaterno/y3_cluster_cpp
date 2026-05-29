// Number counts GPU module WITHOUT selection function.
// Simplified integrand: direct MOR P(lambda|M,z) integration.
//
// This module computes:
//   N = int dlnM int dz  Omega(z) * (dV/dOmega/dz)(z) * n(M,z) * P(lo|M,z) * K(zo|z)
//
// where K(zo|z) is the photo-z kernel (Gaussian smearing).
// Unlike the NumCountsSel module, this does NOT use the pre-tabulated S_ij.
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

#include <iostream>
#include <optional>
#include <vector>

using cosmosis::DataBlock;
using cubacpp::integration_result;

// numberCountsSimple: Direct integration without selection function.
// Integrates over (lo, zt, lnM) with MOR and photo-z kernel.
class numberCountsSimple {
public:
  using grid_t = y3_cluster::grid_t<2>;
  using grid_point_t = grid_t::value_type;

private:
  using volume_t = quad::Volume<double, 3>;

  // Models from datablock
  std::optional<y3_cuda::OMEGA_Z_DES> omega_z_;
  std::optional<y3_cuda::DV_DO_DZ_t> dv_do_dz_;
  std::optional<y3_cuda::HMF_t> hmf_;
  std::optional<y3_cuda::MOR_DES_LOG_t> mor_;
  std::optional<y3_cuda::INT_ZO_ZT_DES_t> int_zo_zt_;

  // Current redshift bin edges
  double zo_low_;
  double zo_high_;

  bool do_cartesian_product_of_bins_ = false;

public:
  explicit numberCountsSimple(cosmosis::DataBlock& cfg)
  {
    auto rc = cfg.get_val(module_label(),
                          "do_cartesian_product_of_bins",
                          false,
                          do_cartesian_product_of_bins_);
    if (rc != DBS_SUCCESS) {
      throw std::runtime_error(
        "numberCountsSimple: failed to read do_cartesian_product_of_bins\n");
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
  }

  void set_grid_point(grid_point_t grid_point)
  {
    zo_low_ = grid_point[0];
    zo_high_ = grid_point[1];
  }

  __host__ __device__ double
  operator()(double lo, double zt, double lnM) const
  {
    // MOR: P(lambda_ob | M, z)
    double const mor_v = (*mor_)(lo, lnM, zt);

    // Common term: Omega(z) * dV/dOmega/dz * dn/dlnM
    double const common = (*omega_z_)(zt) * (*dv_do_dz_)(zt) * (*hmf_)(lnM, zt);

    // Photo-z kernel: integral over z_ob bin
    double const photoz = (*int_zo_zt_)(zo_low_, zo_high_, zt);

    return common * mor_v * photoz;
  }

  static char const* module_label() { return "numberCountsSimple"; }

  static std::vector<volume_t>
  make_integration_volumes(cosmosis::DataBlock& cfg)
  {
    return y3_cuda::make_integration_volumes_wall_of_numbers(
      cfg, numberCountsSimple::module_label(), "lo", "zt", "lnm");
  }

  static grid_t make_grid_points(cosmosis::DataBlock& cfg)
  {
    char const* const label = numberCountsSimple::module_label();
    bool do_cartesian_product_of_bins = false;
    cfg.get_val(label, "do_cartesian_product_of_bins", false,
                do_cartesian_product_of_bins);

    if (do_cartesian_product_of_bins) {
      auto const bottoms = get_vector_double(cfg, label, "zo_bottom");
      auto const widths = get_vector_double(cfg, label, "zo_bin_width");
      if (bottoms.size() != widths.size()) {
        throw std::runtime_error(
          "numberCountsSimple: zo_bottom and zo_bin_width must be same length\n");
      }

      grid_t result;
      result.set_names(std::vector<std::string>{"zo_bottom", "zo_bin_width"});
      for (std::size_t i = 0; i != bottoms.size(); ++i) {
        result.push_back({bottoms[i], bottoms[i] + widths[i]});
      }
      return result;
    }

    return y3_cluster::make_grid_points_wall_of_numbers(
      cfg, numberCountsSimple::module_label(), "zo_low", "zo_high");
  }
};

DEFINE_COSMOSIS_CUDA_INTEGRATION_MODULE(numberCountsSimple)
