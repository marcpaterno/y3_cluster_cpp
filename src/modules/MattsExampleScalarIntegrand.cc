#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "cubacpp/integration_result.hh"
#include "cubacpp/integration_volume.hh"
#include "utils/make_grid_points.hh"
#include "utils/make_integration_volumes.hh"
#include "utils/module_macros.hh"

#include <optional>
#include <vector>

#include "models/hmf_t.hh"


// We write using declarations so that we don't have to type the namespace name
// each time we use these names
using cosmosis::DataBlock;
using cosmosis::ndarray;
using cubacpp::integration_result;

// ExampleScalarIntegrand is a class that models the concept of
// "CosmoSISScalarIntegrand", and is thus suitable for use as the template
// parameter for the class template CosmoSISScalarIntegrationModule.
//
// Notes:
//    1) std::optional<T> is used for data members that are not
//    constructable from CosmoSIS configuration parameters.
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
class ExampleScalarIntegrand {
public:
  // Define the data-type describing a grid point; this should be an
  // instance of std::array<double, N> with N set to the number
  // of different parameters being varied in the grid.
  // The alias we define must be grid_point_t.
  using grid_t = y3_cluster::grid_t<1>;
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
  double sigma_8_;
  std::optional<HMF_t> hmf_;

  // State set for current 'bin' to be integrated.
  double radius_;

public:
  // Initialize my integrand object from the parameters read
  // from the relevant block in the CosmoSIS ini file.
  explicit ExampleScalarIntegrand(cosmosis::DataBlock& config);

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
  double operator()(double x, double y) const;

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


ExampleScalarIntegrand::ExampleScalarIntegrand(DataBlock&)
  : sigma_8_(), hmf_(), radius_()
{}

void
ExampleScalarIntegrand::set_sample(DataBlock& sample)
{
  sigma_8_ = sample.view<double>("cosmological_parameters", "sigma_8");
  // If we had a data member of type std::optional<X>, we would set the
  // value using std::optional::emplace(...) here. emplace takes a set
  // of arguments that it passes to the constructor of X.
}

void
ExampleScalarIntegrand::set_grid_point(grid_point_t const& grid_point)
{
  radius_ = grid_point[0];
}

// This math is totally non-physical and stupid, but it uses all the values
// provided.
double
ExampleScalarIntegrand::operator()(double x, double y) const
{
  // For any data members of type std::optional<X>, we have to use operator*
  // to access the X object (as if we were dereferencing a pointer).
  auto const delta = y - sigma_8_;
  return (x / radius_) + (delta * delta);
}

char const*
ExampleScalarIntegrand::module_label()
{
  return "example_scalar_integrand";
}

// The implementation of make_integration_volumes can be almost the same for
// any CosmoSISIntegrand-type class. Only the names and number of the parameters
// provided need to be changed. It is critical that the names be given in the
// order that correspond to the order of arguments in the class's function call
// operator. While the compiler can verify the number of arguments provided is
// correct, it can not verify that their order matches the order of arguments in
// the function call operator.
std::vector<ExampleScalarIntegrand::volume_t>
ExampleScalarIntegrand::make_integration_volumes(cosmosis::DataBlock& cfg)
{
  return y3_cluster::make_integration_volumes_wall_of_numbers(
    cfg, ExampleScalarIntegrand::module_label(), "x", "y");
}

ExampleScalarIntegrant::grid_t
ExampleScalarIntegrand::make_grid_points(cosmosis::DataBlock& cfg)
{
  return y3_cluster::make_grid_points_cartesian_product(
    cfg, ExampleScalarIntegrand::module_label(), "radii");
}

DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE(ExampleScalarIntegrand)
