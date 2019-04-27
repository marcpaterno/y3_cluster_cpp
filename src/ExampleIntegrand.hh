#ifndef Y3_CLUSTER_EXAMPLE_INTEGRAND_H
#define Y3_CLUSTER_EXAMPLE_INTEGRAND_H

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "cubacpp/integration_volume.hh"
#include "cubacpp/integration_result.hh"

#include <optional>
#include <vector>

// ExampleIntegrand is a class that models the concept of "CosmoSISIntegrand",
// and is thus suitable for use as the template parameter for the class template
// CosmosisIntegrationModule.
//
// Notes:
//    1) std::optional<T> is used for data members that are not
//    constructible from CosmoSIS configuration parameters.
//
//    2) The object as created by the only constructor does not need to be
//    in a callable state.
//
//    3) After a call to set_sample has been made, the object must be in a
//    callable state.
//
//    4) State that *can* be correctly set by the constructor *should* be set by
//    the constructor. Do not needlessly repeat initialization in calls to
//    set_sample.
//
//
class ExampleIntegrand {
private:
  // We define the type alias volume_t to be the right dimensionality
  // of integration volume for our integrand. If we were to change the
  // number of arguments required by the function call operator (below),
  // we would need to also modify this type alias to keep consistent.
  using volume_t = cubacpp::IntegrationVolume<2>;

  // State obtained from configuration. These things should be set in the
  // constructor. Our integrand is a vector-valued function; it will return
  // a vector the same length as radii_, which is the result of evaluating
  // the function to be integrated at that radius.
  std::vector<double> radii_;

  // State obtained from each sample.
  // If there were a type X that did not have a default constructor,
  // we would use std::optional<X> as our data member.
  double sigma_8_;

public:
  // Initialize my integrand object from the parameters read
  // from the relevant block in the CosmoSIS ini file.
  explicit ExampleIntegrand(cosmosis::DataBlock& config);

  // Set any data members from values read from the current sample.
  // Do not attempt to copy the sample!.
  void set_sample(cosmosis::DataBlock& sample);

  // The function to be integrated. All arguments to this function must be of
  // type double, and there must be at least two of them (because our
  // integration routine does not work for functions of one variable). The
  // function is const because calling it does not change the state of the
  // object.
  std::vector<double> operator()(double x, double y) const;

  // finalize_sample() is where products can be put into the cosmosis::DataBlock
  // representing the current sample. The object 'sample' passed to this function
  // will be the exact same object as was passed to the most recent call to
  // set_sample(). The object 'results' will be the results of the integration
  // that has just been done for that sample. This is generally the item which
  // should be put into the sample.
  void finalize_sample(cosmosis::DataBlock& sample,
                       std::vector<cubacpp::integration_results_v> const& results) const;

  // The following non-member (static) function creates a vector of integration
  // volumes (the type alias defined above) based on the parameters read from
  // the configuration block for the module.
  static std::vector<volume_t> make_integration_volumes(
    cosmosis::DataBlock& cfg);
};

#endif
