#include "utils/make_grid_points.hh"
#include "utils/make_integration_volumes.hh"
#include "utils/module_macros.hh"

#include "cosmosis/datablock/datablock.hh"
#include "cubacpp/integration_result.hh"
#include "cubacpp/integration_volume.hh"

#include "models/dv_do_dz_t.hh"
#include "models/hmf_t.hh"
#include "models/int_lc_lt_des_t.hh"
#include "models/int_zo_zt_des_t.hh"
#include "models/lo_lc_t.hh"
#include "models/mor_t.hh"
#include "models/omega_z_des.hh"
#include "models/roffset_t.hh"
#include "models/sig_sum.hh"
#include <optional>
#include <vector>
#include <iostream>
using namespace y3_cluster;
using cosmosis::DataBlock;
using cosmosis::ndarray;
using cubacpp::integration_result;

// MassMiscentY1MortScalarIntegrand is a class that models the concept of
// "CosmoSISScalarIntegrand", and is thus suitable for use as the template
// parameter for the class template CosmoSISScalarIntegrationModule.
//
// Notes:
//    1) std::optional<T> is used for data members that are not
//    constructible from CosmoSIS configuration parameters.
//
//    2) The object as created by the only constructor does not need to be
//    in a callable state.
//
//    3) After calls to both set_sample and set_grid_point have been made, the
//    object must be in a callable state.
//
//    4) State that *can* be correctly set by the constructor *should* be set by
//    the constructor. Do not needlessly repeat initialization in calls to
//    set_sample.
//
//
class MassMiscentY1MortScalarIntegrand {
public:
  // Define the data-type describing a grid point; this should be an
  // instance of std::array<double, N> with N set to the number
  // of different paramaters being varied in the grid.
  // The alias we define must be grid_point_t.
  using grid_point_t = std::array<double, 3>; // we only vary radius.

private:
  // We define the type alias volume_t to be the right dimensionality
  // of integration volume for our integrand. If we were to change the
  // number of arguments required by the function call operator (below),
  // we would need to also modify this type alias to keep consistent.
  using volume_t = cubacpp::IntegrationVolume<4>;

  // State obtained from configuration. These things should be set in the
  // constructor.
  // <none in this example>

  // State obtained from each sample.
  // If there were a type X that did not have a default constructor,
  // we would use std::optional<X> as our data member.
  std::optional<INT_LC_LT_DES_t> lc_lt;
  std::optional<MOR_t> mor;
  std::optional<OMEGA_Z_DES> omega_z;
  std::optional<DV_DO_DZ_t> dv_do_dz;
  std::optional<HMF_t> hmf;
  std::optional<INT_ZO_ZT_DES_t> int_zo_zt;
  std::optional<ROFFSET_t> roffset;
  std::optional<LO_LC_t> lo_lc;

  // State set for current 'bin' to be integrated.
  double zo_low_;
  double zo_high_;
  double rmis_;


public:
  // Initialize my integrand object from the parameters read
  // from the relevant block in the CosmoSIS ini file.
  explicit MassMiscentY1MortScalarIntegrand(cosmosis::DataBlock& config);

  // Set any data members from values read from the current sample.
  // Do not attempt to copy the sample!.
  void set_sample(cosmosis::DataBlock& sample);

  // Set the data for the current bin.
  void set_grid_point(grid_point_t const& pt);

  // The function to be integrated. All arguments to this function must be of
  // type double, and there must be at least two of them (because our
  // integration routine does not work for functions of one variable). The
  // function is const because calling it does not change the state of the
  // object.
  double operator()(
                    double lo,
                    double lc,
                    double zt,
                    double lnM) const;

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
using cubacpp::integration_result;

MassMiscentY1MortScalarIntegrand::MassMiscentY1MortScalarIntegrand(DataBlock&)
  : lc_lt()
  , mor()
  , omega_z()
  , dv_do_dz()
  , hmf()
  , int_zo_zt()
  , roffset()
  , lo_lc()
  , zo_low_()
  , zo_high_()
  , rmis_()
{}

void
MassMiscentY1MortScalarIntegrand::set_sample(DataBlock& sample)
{
  // If we had a data member of type std::optional<X>, we would set the
  // value using std::optional::emplace(...) here. emplace takes a set
  // of arguments that it passes to the constructor of X.
  lc_lt.emplace(sample);
  mor.emplace(sample);
  dv_do_dz.emplace(sample);
  hmf.emplace(sample);
  omega_z.emplace(sample);
  roffset.emplace(sample);
  lo_lc.emplace(sample);
}

void
MassMiscentY1MortScalarIntegrand::set_grid_point(grid_point_t const& grid_point)
{
  zo_low_ = grid_point[0];
  zo_high_ = grid_point[1];
  rmis_ = grid_point[2];
}

double
MassMiscentY1MortScalarIntegrand::operator()(double lo,
                             double lc,
                             double zt,
                             double lnM) const
{
  // For any data members of type std::optional<X>, we have to use operator*
  // to access the X object (as if we were dereferencing a pointer).
  double mass = std::exp(lnM);
  double common_term = (*lo_lc)(lo, lc, rmis_) *
                       (*mor)(lc, lnM, zt) *
                       (*dv_do_dz)(zt) * (*hmf)(lnM, zt) * (*omega_z)(zt);
  auto const val = mass*(*int_zo_zt)(zo_low_, zo_high_, zt) * common_term;
  return val;
}

char const*
MassMiscentY1MortScalarIntegrand::module_label()
{
  return "MassMiscentY1MortScalarIntegrand";
}

// The implementation of make_integration_volumes can be almost the same for
// any CosmoSISIntegrand-type class. Only the names and number of the parameters
// provided need to be changed. It is critical that the names be given in the
// order that correspond to the order of arguments in the class's function call
// operator. While the compiler can verify the number of arguments provided is
// correct, it can not verify that their order matches the order of arguments in
// the function call operator.
std::vector<MassMiscentY1MortScalarIntegrand::volume_t>
MassMiscentY1MortScalarIntegrand::make_integration_volumes(cosmosis::DataBlock& cfg)
{
  return y3_cluster::make_integration_volumes_wall_of_numbers(
    cfg, MassMiscentY1MortScalarIntegrand::module_label(), "lo", "lc", "zt", "lnm");
}

std::vector<MassMiscentY1MortScalarIntegrand::grid_point_t>
MassMiscentY1MortScalarIntegrand::make_grid_points(cosmosis::DataBlock& cfg)
{
  return y3_cluster::make_grid_points_cartesian_product(
    cfg, MassMiscentY1MortScalarIntegrand::module_label(), "zo_low", "zo_high", "rmis_array");
}

DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE(MassMiscentY1MortScalarIntegrand)
