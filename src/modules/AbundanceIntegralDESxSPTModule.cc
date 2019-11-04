#include "cosmosis/datablock/ndarray.hh"
#include "utils/make_grid_points.hh"
#include "utils/make_integration_volumes.hh"
#include "utils/module_macros.hh"

#include "cosmosis/datablock/datablock.hh"
#include "cubacpp/integration_result.hh"
#include "cubacpp/integration_volume.hh"

#include "models/zo_zt_des_t.hh"
#include "models/sptxdes/omega_z_y3xspt.hh"
#include "models/dv_do_dz_t.hh"
#include "models/hmf_t.hh"
#include "models/lc_lt_y1_t.hh"
#include "models/sptxdes/pxizeta_t.hh"
#include "models/sptxdes/mor_y3xspt_t.hh"

#include <optional>
#include <vector>
#include <cmath>

using namespace y3_cluster;

class AbundanceIntegralDesXSpt{
public:
  // 3 varied grid parameters: lo, xi, zo
  using grid_point_t = std::array<double, 4>;

private:
  // Number of integration dimensions
  using volume_t = cubacpp::IntegrationVolume<4>;

  // State obtained from configuration
  ZO_ZT_DES_t p_zo_zt;
  OMEGA_Z_Y3XSPT omega_z;
  LC_LT_Y1_t p_lo_ltzt;
  PXiZeta_t p_xi_zeta;

  // State obtained from each sample.
  // If there were a type X that did not have a default constructor,
  // we would use std::optional<X> as our data member.
  std::optional<DV_DO_DZ_t> dv_dodz;
  std::optional<HMF_t> hmf;
  std::optional<MOR_Y3xSPT_t> p_ltzeta_lnmzt;

  // State set for current 'bin' to be integrated.
  double lamobs_;
  double xi_;
  double zobs_;
  double gamma_field_;

public:
  // Initialize my integrand object from the parameters read
  // from the relevant block in the CosmoSIS ini file.
  explicit AbundanceIntegralDesXSpt(cosmosis::DataBlock& config);

  // Set any data members from values read from the current sample.
  // Do not attempt to copy the sample!.
  void set_sample(cosmosis::DataBlock& sample);

  // Set the data for the current bin.
  void set_grid_point(grid_point_t const& grid_point);

  // The function to be integrated. All arguments to this function must be of
  // type double, and there must be at least two of them (because our
  // integration routine does not work for functions of one variable). The
  // function is const because calling it does not change the state of the
  // object.
  double operator()(double ztrue, double lnM200m, double lamtrue, double zeta) const;

  // module_label() is a non-member (static) function that returns the label for
  // this module. The name this returns
  // is the name that must be used in the 'ini file' for configuring the module
  // made with this class.
  // We return char const* rather than std::string to avoid some needless memory
  // allocations.
  static char const* module_label();

  // The following non-member (static) function creates a vector of integration
  // volumes (the type alias defined above) based on the parameters read from
  // the configuration block for the module.
  static std::vector<volume_t> make_integration_volumes(
    cosmosis::DataBlock& cfg);

  // The following non-member (static) function creates a vector of grid points
  // on which the integration results are to be evaluated, based on parameters
  // read from the configuration block for the module.
  static std::vector<grid_point_t> make_grid_points(cosmosis::DataBlock& cfg);
};


// We write using declarations so that we don't have to type the namespace name
// each time we use these names
using cosmosis::DataBlock;
using cosmosis::ndarray;
using cubacpp::integration_result;

// Set up integral that does not depend on sample
AbundanceIntegralDesXSpt::AbundanceIntegralDesXSpt(DataBlock&)
  : p_zo_zt()
  , omega_z()
  , p_lo_ltzt()
  , p_xi_zeta()
  , dv_dodz()
  , hmf()
  , p_ltzeta_lnmzt()
  , lamobs_()
  , xi_()
  , zobs_()
  , gamma_field_()
{}

// Set up the integration using the current sample
void
AbundanceIntegralDesXSpt::set_sample(DataBlock& sample)
{
  dv_dodz.emplace(sample);
  hmf.emplace(sample);
  p_ltzeta_lnmzt.emplace(sample);
}

// Set the grid points for the current integration
void
AbundanceIntegralDesXSpt::set_grid_point(grid_point_t const& grid_point)
{
  lamobs_ = grid_point[0];
  xi_ = grid_point[1];
  zobs_ = grid_point[2];
  gamma_field_ = grid_point[3];
}

// Evaluate the integrand
// Anything initilized with std::optional needs an asterisk
double
AbundanceIntegralDesXSpt::operator()(double ztrue,
                                     double lnM200m,
                                     double lamtrue,
                                     double zeta) const
{
  return p_zo_zt(zobs_, ztrue)
         * omega_z(ztrue)
         * (*dv_dodz)(ztrue)
         * (*hmf)(lnM200m, ztrue)
         * p_lo_ltzt(lamobs_, lamtrue, ztrue)
         * p_xi_zeta(xi_, zeta)
         * (*p_ltzeta_lnmzt)(lamtrue, zeta, ztrue, lnM200m, gamma_field_);
}

// Name of the module in the datablock
char const*
AbundanceIntegralDesXSpt::module_label()
{
  return "abundance_integral";
}

// Set the bounds of integration
std::vector<AbundanceIntegralDesXSpt::volume_t>
AbundanceIntegralDesXSpt::make_integration_volumes(cosmosis::DataBlock& cfg)
{
  return y3_cluster::make_integration_volumes_wall_of_numbers(
    cfg, AbundanceIntegralDesXSpt::module_label(), "ztrue",
                                                   "lnM200m",
                                                   "lamtrue",
                                                   "zeta");
}

// Set the grid points to evaluate the integral at
std::vector<AbundanceIntegralDesXSpt::grid_point_t>
AbundanceIntegralDesXSpt::make_grid_points(cosmosis::DataBlock& cfg)
{
  return y3_cluster::make_grid_points_wall_of_numbers(
    cfg, AbundanceIntegralDesXSpt::module_label(), "lamobs",
                                                   "xi",
                                                   "zobs",
                                                   "gamma_field");
}


DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE(AbundanceIntegralDesXSpt);
