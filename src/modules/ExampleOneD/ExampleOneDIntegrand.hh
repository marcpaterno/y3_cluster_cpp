#ifndef Y3_CLUSTER_EXAMPLE_ONED_INTEGRAND_H
#define Y3_CLUSTER_EXAMPLE_ONED_INTEGRAND_H

#include "cosmosis/datablock/datablock.hh"

#include "utils/make_grid_points.hh"

#include <optional>
#include <vector>
#include <string>

// ExampleOneDIntegrand is a class that models the concept of
// "OneDIntegration", and is thus suitable for use as the template
// parameter for the class template OneDIntegrationModule.
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

// In this example, we will build a class to integrate f(x) = a*x + b + c
// a will be on a grid loaded from a configuration file
// b will be fixed but loaded from a configuration file
// c will be loaded from the sample

class ExampleOneDIntegrand {
public:
  using grid_t = y3_cluster::grid_t<1>;
  using grid_point_t = grid_t::value_type;

private:
  // State obtained from configuration. Set in constructor
  double b_;

  // State obtained from each sample. Set in set_sample
  // If there were a type X that did not have a default constructor,
  // we would use std::optional<X> as our data member.
  double c_;

  // State obtained from grid. Set in set_grid_point
  double a_;

public:
  // module_label() is a non-member (static) function that returns the label for
  // this module. The name this returns
  // is the name that must be used in the 'ini file' for configuring the module
  // made with this class.
  // We return char const* rather than std::string to avoid some needless memory
  // allocations.
  static char const* module_label();

  // a non-member function that returns a string containing the variable of
  // integration to reference in the configuration file
  static char const* integration_variable();

  // The function to be integrated. This function should take 
  //
  //
  // All arguments to this function must be of
  // type double, and there must be at least two of them (because our
  // integration routine does not work for functions of one variable). The
  // function is const because calling it does not change the state of the
  // object.
  double operator()(double x) const;

  // Initialize my integrand object from the parameters read
  // from the relevant block in the CosmoSIS ini file.
  explicit ExampleOneDIntegrand(cosmosis::DataBlock& config);

  // --------------------------------------------------------------------------
  // --------- Prototypes below will be unchanged for most users --------------
  // --------------------------------------------------------------------------
  // Set any data members from values read from the current sample.
  // Do not attempt to copy the sample!.
  void set_sample(cosmosis::DataBlock& sample);

  // Set the data for the current bin.
  void set_grid_point(grid_point_t const& pt);

  // The following non-member (static) function creates a vector of integration
  // volumes (the type alias defined above) based on the parameters read from
  // the configuration block for the module.
  static std::vector<std::array<double, 2>> make_integration_volumes(
    cosmosis::DataBlock& cfg);

  // The following non-member (static) function creates a vector of grid points
  // on which the integration results are to be evaluated, based on parameters
  // read from the configuration block for the module.
  static grid_t make_grid_points(cosmosis::DataBlock& cfg);
};

#endif
