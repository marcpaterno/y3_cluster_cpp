#include "cosmosis/datablock/ndarray.hh"
#include "utils/make_grid_points.hh"
#include "utils/make_integration_volumes.hh"
#include "utils/module_macros.hh"
#include "utils/datablock_reader.hh"

#include "cosmosis/datablock/datablock.hh"
#include "cubacpp/integration_result.hh"
#include "cubacpp/integration_volume.hh"

#include "models/int_zo_zt_des_t.hh"
#include "models/sptxdes/omega_z_y3xspt.hh"
#include "models/dv_do_dz_t.hh"
#include "models/hmf_t.hh"
#include "models/lc_lt_y1_t.hh"
#include "models/sptxdes/int_pxizeta.hh"
#include "models/sptxdes/mor_y3xspt_t.hh"

#include <optional>
#include <vector>
#include <string>
#include <cmath>

using namespace y3_cluster;

class SelectionIntegralDESxSPT{
public:
  // Grid parameters
  using grid_point_t = std::array<double, 2>;

private:
  // Number of integration dimensions
  using volume_t = cubacpp::IntegrationVolume<5>;

  // State obtained from configuration
  INT_ZO_ZT_DES_t int_p_zo_zt;
  OMEGA_Z_Y3XSPT omega_z;
  LC_LT_Y1_t p_lo_ltzt;
  INT_PXIZETA int_pxizeta;

  // State obtained from each sample.
  // If there were a type X that did not have a default constructor,
  // we would use std::optional<X> as our data member.
  std::optional<DV_DO_DZ_t> dv_dodz;
  std::optional<HMF_t> hmf;
  std::optional<MOR_Y3xSPT_t> p_ltzeta_lnmzt;

  // State set for current 'bin' to be integrated.
  double zobs_min_;
  double zobs_max_;

public:
  // Initialize my integrand object from the parameters read
  // from the relevant block in the CosmoSIS ini file.
  explicit SelectionIntegralDESxSPT(cosmosis::DataBlock& config);

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
  double operator()(double lamobs, double lamtrue, double lnM200m,
                    double ztrue, double zeta) const;

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
SelectionIntegralDESxSPT::SelectionIntegralDESxSPT(DataBlock&)
  : int_p_zo_zt()
  , omega_z()
  , p_lo_ltzt()
  , int_pxizeta()
  , dv_dodz()
  , hmf()
  , p_ltzeta_lnmzt()
  , zobs_min_()
  , zobs_max_()
{}

// Set up the integration using the current sample
void
SelectionIntegralDESxSPT::set_sample(DataBlock& sample)
{
  dv_dodz.emplace(sample);
  hmf.emplace(sample);
  p_ltzeta_lnmzt.emplace(sample);
}

// Set the grid points for the current integration
void
SelectionIntegralDESxSPT::set_grid_point(grid_point_t const& grid_point)
{
  zobs_min_ = grid_point[0];
  zobs_max_ = grid_point[1];
}

// Evaluate the integrand
// Anything initilized with std::optional needs an asterisk
double
SelectionIntegralDESxSPT::operator()(double lamobs,
                                     double lamtrue,
                                     double lnM200m,
                                     double ztrue,
                                     double zeta) const
{
  double gamma_field = 1.0;
  return int_p_zo_zt(zobs_min_, zobs_max_, ztrue)
         * omega_z(ztrue)
         * (*dv_dodz)(ztrue)
         * (*hmf)(lnM200m, ztrue)
         * p_lo_ltzt(lamobs, lamtrue, ztrue)
         * int_pxizeta(zeta)
         * (*p_ltzeta_lnmzt)(lamtrue, zeta, ztrue, lnM200m, gamma_field);
}

// Name of the module in the datablock
char const*
SelectionIntegralDESxSPT::module_label()
{
  return "selection_integral";
}


// Set the bounds of integration
std::vector<SelectionIntegralDESxSPT::volume_t>
SelectionIntegralDESxSPT::make_integration_volumes(cosmosis::DataBlock& cfg)
{
  return y3_cluster::make_integration_volumes_wall_of_numbers(
      cfg, SelectionIntegralDESxSPT::module_label(), "lamobs",
                                                     "lamtrue",
                                                     "lnM200m",
                                                     "ztrue",
                                                     "zeta");
}

// Set the grid points to evaluate the integral at
std::vector<SelectionIntegralDESxSPT::grid_point_t>
SelectionIntegralDESxSPT::make_grid_points(cosmosis::DataBlock& cfg)
{
  return y3_cluster::make_grid_points_wall_of_numbers(
    cfg, SelectionIntegralDESxSPT::module_label(), "zobs_min",
                                                   "zobs_max");
}


DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE(SelectionIntegralDESxSPT)
