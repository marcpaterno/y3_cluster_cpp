#include "utils/make_grid_points.hh"
#include "utils/cuda_module_macros.cuh"

#include "cosmosis/datablock/datablock.hh"
#include "cubacpp/integration_result.hh"
#include "cudaPagani/quad/util/Volume.cuh"

#include "models/dv_do_dz_t.cuh"
#include "models/lo_lc_t.cuh"
#include "models/roffset_t.cuh"

#include <iostream>
#include <optional>
#include <vector>

using namespace y3_cluster;
using cosmosis::DataBlock;
using cosmosis::ndarray;
using cubacpp::integration_result;

// This is a class that models the concept of
// "CosmoSISScalarIntegrand", and is thus suitable for use as the template
// parameter for the class template CosmoSISSICUDAModule.
//
// CAUTION: This is *very* preliminary. Don't hurt yourself.
//
class SigmaMiscentY1ScalarIntegrand {
public:
  using grid_t = y3_cluster::grid_t<1>;
  using grid_point_t = grid_t::value_type;

private:
  // We define the type alias volume_t to be the right dimensionality
  // of integration volume for our integrand. If we were to change the
  // number of arguments required by the function call operator (below),
  // we would need to also modify this type alias to keep consistent.
  using volume_t = quad::Volume<double, 3>;

  // State obtained from configuration. These things should be set in the
  // constructor.
  // <none in this example>

  // State obtained from each sample.
  // If there were a type X that did not have a default constructor,
  // we would use std::optional<X> as our data member.
  std::optional<DV_DO_DZ_t> dv_do_dz;
  std::optional<LO_LC_t> lo_lc;
  std::optional<ROFFSET_t> roffset;

  // State set for current 'bin' to be integrated.
  double radius_;

public:
  // Initialize my integrand object from the parameters read
  // from the relevant block in the CosmoSIS ini file.
  explicit SigmaMiscentY1ScalarIntegrand(cosmosis::DataBlock& config);

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
  __device__
  double operator()(double lo,
                    double lc,
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

SigmaMiscentY1ScalarIntegrand::SigmaMiscentY1ScalarIntegrand(DataBlock&)
  : dv_do_dz()
  , roffset()
  , lo_lc()
  , radius_()
{}

void
SigmaMiscentY1ScalarIntegrand::set_sample(DataBlock& sample)
{
  // If we had a data member of type std::optional<X>, we would set the
  // value using std::optional::emplace(...) here. emplace takes a set
  // of arguments that it passes to the constructor of X.
  dv_do_dz.emplace(sample);
  roffset.emplace(sample);
  lo_lc.emplace(sample);
}

void
SigmaMiscentY1ScalarIntegrand::set_grid_point(grid_point_t const& grid_point)
{
  radius_ = grid_point[0];
}

__device__ double
SigmaMiscentY1ScalarIntegrand::operator()(double lo,
                                          double lc,
                                          double rmis,
                                          double zt) const
{
  double common_term = (*roffset)(rmis) * (*lo_lc)(lo, lc, rmis) *
                       (*dv_do_dz)(zt) / 2.0 / 3.1415926535897;
  double scaled_Rmis = sqrt(radius_ * radius_ + rmis * rmis +
                                 2 * rmis * radius_;
  auto const val = scaled_Rmis * common_term;
  return val;
}

char const*
SigmaMiscentY1ScalarIntegrand::module_label()
{
  return "SigmaMiscentY1ScalarIntegrand";
}

// The implementation of make_integration_volumes can be almost the same for
// any CosmoSISIntegrand-type class. Only the names and number of the parameters
// provided need to be changed. It is critical that the names be given in the
// order that correspond to the order of arguments in the class's function call
// operator. While the compiler can verify the number of arguments provided is
// correct, it can not verify that their order matches the order of arguments in
// the function call operator.
std::vector<SigmaMiscentY1ScalarIntegrand::volume_t>
SigmaMiscentY1ScalarIntegrand::make_integration_volumes(
  cosmosis::DataBlock& cfg)
{
  return y3_cluster::make_integration_volumes_wall_of_numbers(
    cfg,
    SigmaMiscentY1ScalarIntegrand::module_label(),
    "lo",
    "lc",
    "lt",
    "zt",
    "lnm",
    "rmis",
    "theta");
}

SigmaMiscentY1ScalarIntegrand::grid_t
SigmaMiscentY1ScalarIntegrand::make_grid_points(cosmosis::DataBlock& cfg)
{
  return y3_cluster::make_grid_points_cartesian_product(
    cfg,
    SigmaMiscentY1ScalarIntegrand::module_label(),
    "zo_low",
    "zo_high",
    "radii");
}

DEFINE_COSMOSIS_CUDA_INTEGRATION_MODULE(SigmaMiscentY1ScalarIntegrand)
