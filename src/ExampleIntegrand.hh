#ifndef Y3_CLUSTER_EXAMPLE_INTEGRAND_H
#define Y3_CLUSTER_EXAMPLE_INTEGRAND_H

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "hmf_t.hh"
#include "integration_range.hh"

#include <optional>

// ExampleIntegrand is a class that models the concept of "CosmoSISIntegrand",
// and is suitable for use as the template parameter for the class template
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
  // State obtained from configuration. These things should be set in the
  // constructor.
  std::vector<double> Rs_;
  std::vector<y3_cluster::IntegrationRange> lnM_ranges_;
  std::vector<y3_cluster::IntegrationRange> z_ranges_;

  // State obtained from each sample. We use std::optional to hold the function
  // object(s), so that we do not need to require then to be
  // default-constructible and assignable.
  double sigma_8_;
  std::optional<y3_cluster::HMF_t> hmf_;

public:
  // Initialize my integrand object from the parameters read
  // from the relevant block in the CosmoSIS ini file.
  explicit ExampleIntegrand(cosmosis::DataBlock& config);

  // Set any data members from values read from the current sample.
  void set_sample(cosmosis::DataBlock& sample);

  // The function to be integrated. All arguments to this function must be of
  // type double, and there must be at least two of them (because our
  // integration routine does not work for functions of one variable). The
  // function is const because calling it does not change the state of the
  // object.
  std::vector<double> execute(double z, double lambda) const;
};

#endif