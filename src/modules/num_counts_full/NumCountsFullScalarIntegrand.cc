// Full-integral reference for number counts N_i[1] per richness bin.
//
//   N_i[1] = ∫∫∫ dlt dz dlnM  n(M,z) · dV/dΩdz · Ω(z)
//                              · B_i(lt, z) · P(lt | M, z; φ)
//
// Brute-force reference for N_i[1]: integrates the full triple integral
// directly so downstream tests can validate any approximation that lands
// on the same datablock layout.
//
// Scope is intentionally narrower than the Y1 forward model — no
// miscentering, no photo-z.
#include "utils/make_grid_points.hh"
#include "utils/make_integration_volumes.hh"
#include "utils/module_macros.hh"

#include "cosmosis/datablock/datablock.hh"
#include "cubacpp/integration_volume.hh"

#include "models/dv_do_dz_t.hh"
#include "models/hmf_t.hh"
#include "models/mor_hod_t.hh"
#include "models/omega_z_des.hh"
#include "models/projection_y3_b_i_t.hh"

#include <array>
#include <optional>
#include <vector>

class NumCountsFullScalarIntegrand {
public:
  using grid_t = y3_cluster::grid_t<1>;
  using grid_point_t = grid_t::value_type;

private:
  using volume_t = cubacpp::IntegrationVolume<3>;

  std::optional<y3_cluster::HMF_t> hmf_;
  std::optional<y3_cluster::DV_DO_DZ_t> dv_do_dz_;
  std::optional<y3_cluster::OMEGA_Z_DES> omega_z_;
  std::optional<y3_cluster::MOR_HOD_t> mor_;
  std::array<std::optional<y3_cluster::PROJ_Y3_B_i_t>, 4> b_i_;

  int current_bin_{0};

public:
  explicit NumCountsFullScalarIntegrand(cosmosis::DataBlock& /*cfg*/)
  {
    for (int i = 0; i < 4; ++i) { b_i_[i].emplace(i); }
  }

  void
  set_sample(cosmosis::DataBlock& sample)
  {
    hmf_.emplace(sample);
    dv_do_dz_.emplace(sample);
    omega_z_.emplace(sample);
    mor_.emplace(sample);
  }

  void
  set_grid_point(grid_point_t const& pt)
  {
    current_bin_ = static_cast<int>(pt[0]);
  }

  double
  operator()(double lt, double zt, double lnM) const
  {
    return (*hmf_)(lnM, zt) * (*dv_do_dz_)(zt) * (*omega_z_)(zt) *
           (*b_i_[current_bin_])(lt, zt) * (*mor_)(lt, lnM, zt);
  }

  static char const*
  module_label()
  {
    return "NumCountsFullScalarIntegrand";
  }

  static std::vector<volume_t>
  make_integration_volumes(cosmosis::DataBlock& cfg)
  {
    return y3_cluster::make_integration_volumes_wall_of_numbers(
      cfg, module_label(), "lt", "zt", "lnm");
  }

  static grid_t
  make_grid_points(cosmosis::DataBlock& cfg)
  {
    return y3_cluster::make_grid_points_wall_of_numbers(
      cfg, module_label(), "lambda_bin");
  }
};

DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE(NumCountsFullScalarIntegrand)
