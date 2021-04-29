#include "utils/make_grid_points.hh"
#include "utils/make_integration_volumes.hh"
#include "utils/module_macros.hh"

#include "cosmosis/datablock/datablock.hh"
#include "cubacpp/integration_result.hh"
#include "cubacpp/integration_volume.hh"

#include "models/hmf_t.hh"
#include "models/sig_sum.hh"

#include <iostream>
#include <optional>
#include <vector>
using namespace y3_cluster;
using cosmosis::DataBlock;
using cosmosis::ndarray;
using cubacpp::integration_result;

// SnapshotScalarNCIntegrand is a class that models the concept of
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
class SnapshotScalarNCIntegrand {
public:
  using grid_t = y3_cluster::grid_t<2>;
  using grid_point_t = grid_t::value_type;

private:
  // We define the type alias volume_t to be the right dimensionality
  // of integration volume for our integrand. If we were to change the
  // number of arguments required by the function call operator (below),
  // we would need to also modify this type alias to keep consistent.
  using volume_t = cubacpp::IntegrationVolume<2>;

  // State obtained from configuration. These things should be set in the
  // constructor.
  // <none in this example>

  // State obtained from each sample.
  // If there were a type X that did not have a default constructor,
  // we would use std::optional<X> as our data member.
  std::optional<HMF_t> hmf;
  std::optional<SIG_SUM> sigma;

  // State set for current 'bin' to be integrated.
  double radius_;
  double zt_;

public:
  // Initialize my integrand object from the parameters read
  // from the relevant block in the CosmoSIS ini file.
  explicit SnapshotScalarNCIntegrand(cosmosis::DataBlock& config);

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
  double operator()(double lt, double lnM) const;

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

SnapshotScalarNCIntegrand::SnapshotScalarNCIntegrand(DataBlock&)
  : hmf(), sigma(), radius_(), zt_()
{}

void
SnapshotScalarNCIntegrand::set_sample(DataBlock& sample)
{
  // If we had a data member of type std::optional<X>, we would set the
  // value using std::optional::emplace(...) here. emplace takes a set
  // of arguments that it passes to the constructor of X.
  hmf.emplace(sample);
  sigma.emplace(sample);
}

void
SnapshotScalarNCIntegrand::set_grid_point(grid_point_t const& grid_point)
{
  radius_ = grid_point[1];
  zt_ = grid_point[0];
}

double
SnapshotScalarNCIntegrand::operator()(double /* lt */, double lnM) const
{
  // For any data members of type std::optional<X>, we have to use operator*
  // to access the X object (as if we were dereferencing a pointer).
  auto constexpr simulation_cosmic_volume = 165.0 * 165.0 * 165.0;
  auto const val = simulation_cosmic_volume
                   //* (*mor)(lt, lnM, zt_) * (*hmf)(lnM, zt_);
                   * (*hmf)(lnM, zt_);
  //* (*sigma)(radius_, lnM, zt_);
  return val;
}

char const*
SnapshotScalarNCIntegrand::module_label()
{
  return "SnapshotScalarNCIntegrand";
}

// The implementation of make_integration_volumes can be almost the same for
// any CosmoSISIntegrand-type class. Only the names and number of the parameters
// provided need to be changed. It is critical that the names be given in the
// order that correspond to the order of arguments in the class's function call
// operator. While the compiler can verify the number of arguments provided is
// correct, it can not verify that their order matches the order of arguments in
// the function call operator.
std::vector<SnapshotScalarNCIntegrand::volume_t>
SnapshotScalarNCIntegrand::make_integration_volumes(cosmosis::DataBlock& cfg)
{
  return y3_cluster::make_integration_volumes_wall_of_numbers(
    cfg, SnapshotScalarNCIntegrand::module_label(), "lt", "lnm");
}

SnapshotScalarNCIntegrand::grid_t
SnapshotScalarNCIntegrand::make_grid_points(cosmosis::DataBlock& cfg)
{
  return y3_cluster::make_grid_points_cartesian_product(
    cfg, SnapshotScalarNCIntegrand::module_label(), "snapshot_zs", "radii");
}

DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE(SnapshotScalarNCIntegrand)
