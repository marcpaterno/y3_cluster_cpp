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
#include "models/mor_des_log_t.hh"
#include "models/omega_z_des.hh"
#include "models/roffset_t.hh"
#include "models/kappa_max.hh"

#include <iostream>
#include <optional>
#include <vector>

using cosmosis::DataBlock;
using cosmosis::ndarray;
using cubacpp::integration_result;


// This is a class that models the concept of
// "CosmoSISCUDAScalarIntegrand", and is thus suitable for use as the template
// parameter for the class template CosmoSISScalarIntegrationModule.
//
class avgMiscApprox1CardCPU {
public:
  using grid_t = y3_cluster::grid_t<3>;
  using grid_point_t = grid_t::value_type;

private:
  // We define the type alias volume_t to be the right dimensionality
  // of integration volume for our integrand. If we were to change the
  // number of arguments required by the function call operator (below),
  // we would need to also modify this type alias to keep consistent.
  using volume_t = cubacpp::IntegrationVolume<6>;

  // State obtained from configuration. These things should be set in the
  // constructor.
  // <none in this example>

  // State obtained from each sample.
  // If there were a type X that did not have a default constructor,
  // we would use std::optional<X> as our data member.
  //
  // the volume
  // std::optional<y3_cluster::MOR_t> mor;
  std::optional<y3_cluster::INT_LC_LT_DES_t> lc_lt;
  std::optional<y3_cluster::MOR_DES_LOG_t> mor;
  std::optional<y3_cluster::OMEGA_Z_DES> omega_z;
  std::optional<y3_cluster::DV_DO_DZ_t> dv_do_dz;
  std::optional<y3_cluster::HMF_t> hmf;
  std::optional<y3_cluster::INT_ZO_ZT_DES_t> int_zo_zt;
  std::optional<y3_cluster::ROFFSET_t> roffset;
  std::optional<y3_cluster::KAPPA_MAX> gamma;

  // State set for current 'bin' to be integrated.
  double zo_low_;
  double zo_high_;
  double radius_;

  bool do_cartesian_product_of_bins_;
public:
  // Initialize my integrand object from the parameters read
  // from the relevant block in the CosmoSIS ini file.
  explicit avgMiscApprox1CardCPU(cosmosis::DataBlock& config);

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
  double operator()(double lo,
                    double lt,
                    double zt,
                    double lnM,
                    double rmis,
                    double theta) const;

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
  static grid_t make_grid_points(cosmosis::DataBlock& cfg);
};

// We write using declarations so that we don't have to type the namespace name
// each time we use these names
using cosmosis::DataBlock;
using cubacpp::integration_result;

avgMiscApprox1CardCPU::avgMiscApprox1CardCPU( DataBlock& cfg)
{
  auto rc = 
    cfg.get_val(module_label(),
                "do_cartesian_product_of_bins",
                false,
                do_cartesian_product_of_bins_);
  if (rc != DBS_SUCCESS) {
    throw std::runtime_error("summedNumbersCentY1 failed to find do_cartesian_product_of_bins\n");
  }
}

void
avgMiscApprox1CardCPU::set_sample(DataBlock& sample)
{
  // If we had a data member of type std::optional<X>, we would set the
  // value using std::optional::emplace(...) here. emplace takes a set
  // of arguments that it passes to the constructor of X.
  omega_z.emplace(sample);
  dv_do_dz.emplace(sample);
  hmf.emplace(sample);
  mor.emplace(sample);
  gamma.emplace(sample);
  lc_lt.emplace(sample);
  roffset.emplace(sample);
}

void
avgMiscApprox1CardCPU::set_grid_point(
  grid_point_t const& grid_point)
{
  radius_ = grid_point[2];
  zo_low_ = grid_point[0];
  zo_high_ = grid_point[1];
}

double
avgMiscApprox1CardCPU::operator()(double lo,
                                  double lt,
                                  double zt,
                                  double lnM,
                                  double rmis,
                                  double theta) const
{
  // For any data members of type std::optional<X>, we have to use operator*
  // to access the X object (as if we were dereferencing a pointer).
  
  // Mis-centering approximation; lambda is not affected by rmis
  double const lc = lo;

  double const roffset_v = (*roffset)(rmis);
  double scaled_Rmis = std::sqrt(radius_ * radius_ + rmis * rmis +
                                 2 * rmis * radius_ * std::cos(theta));

  double const mor_v = (*mor)(lo, lnM, zt);
  double common_term = (*omega_z)(zt) * (*dv_do_dz)(zt) * (*hmf)(lnM, zt) * mor_v ;
  auto const val = (*gamma)(scaled_Rmis, lnM, zt) * roffset_v * (*lc_lt)(lc, lt, zt) *
                   (*int_zo_zt)(zo_low_, zo_high_, zt) * common_term;
  return val;
}

// string must match section block in pipeline.ini file
char const*
avgMiscApprox1CardCPU::module_label()
{
  return "gammaMisc";
}

// The implementation of make_integration_volumes can be almost the same for
// any CosmoSISIntegrand-type class. Only the names and number of the parameters
// provided need to be changed. It is critical that the names be given in the
// order that correspond to the order of arguments in the class's function call
// operator. While the compiler can verify the number of arguments provided is
// correct, it can not verify that their order matches the order of arguments in
// the function call operator.
std::vector<avgMiscApprox1CardCPU::volume_t>
avgMiscApprox1CardCPU::make_integration_volumes(
  cosmosis::DataBlock& cfg)
{
  return y3_cluster::make_integration_volumes_wall_of_numbers(
    cfg, avgMiscApprox1CardCPU::module_label(),     
    "lo",
    "lt",
    "zt",
    "lnm",
    "rmis",
    "theta");
}

avgMiscApprox1CardCPU::grid_t
avgMiscApprox1CardCPU::make_grid_points(cosmosis::DataBlock& cfg)
{
  return y3_cluster::make_grid_points_wall_of_numbers(
    cfg,
    avgMiscApprox1CardCPU::module_label(),
    "zo_low",
    "zo_high",
    "radii");
}

DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE(avgMiscApprox1CardCPU)
