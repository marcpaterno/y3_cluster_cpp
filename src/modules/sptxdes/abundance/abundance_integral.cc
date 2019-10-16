#include "abundance_integral.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "utils/make_grid_points.hh"
#include "utils/make_integration_volumes.hh"

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
         * p_xi_zeta(xi_, zeta, gamma_field_)
         * (*p_ltzeta_lnmzt)(lamtrue, zeta, ztrue, lnM200m);
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










