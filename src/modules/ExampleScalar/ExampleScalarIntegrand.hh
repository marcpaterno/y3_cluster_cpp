#ifndef Y3_CLUSTER_EXAMPLE_SCALAR_INTEGRAND_H
#define Y3_CLUSTER_EXAMPLE_SCALAR_INTEGRAND_H

#include "cosmosis/datablock/datablock.hh"
#include "cubacpp/integration_result.hh"
#include "cubacpp/integration_volume.hh"

#include "utils/make_grid_points.hh"

#include <optional>
#include <vector>

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

#endif
