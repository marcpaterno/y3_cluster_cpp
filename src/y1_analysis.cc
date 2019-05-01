#include "make_integration_volumes.hh"

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "cubacpp/integration_volume.hh"
#include "cubacpp/integration_result.hh"

#include <optional>
#include <vector>
#include "lc_lt_t.hh"
using namespace y3_cluster;

// y1_analysis is a class that models the concept of "CosmoSISIntegrand",
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
class y1_analysis {
private:
  // We define the type alias volume_t to be the right dimensionality
  // of integration volume for our integrand. If we were to change the
  // number of arguments required by the function call operator (below),
  // we would need to also modify this type alias to keep consistent.
  using volume_t = cubacpp::IntegrationVolume<3>;

  // State obtained from each sample.
  // If there were a type X that did not have a default constructor,
  // we would use std::optional<X> as our data member.
  std::optional<LC_LT_t> lc_lt;

public:
  // Initialize my integrand object from the parameters read
  // from the relevant block in the CosmoSIS ini file.
  explicit y1_analysis(cosmosis::DataBlock& config);

  // Set any data members from values read from the current sample.
  // Do not attempt to copy the sample!.
  void set_sample(cosmosis::DataBlock& sample);

  // The function to be integrated. All arguments to this function must be of
  // type double, and there must be at least two of them (because our
  // integration routine does not work for functions of one variable). The
  // function is const because calling it does not change the state of the
  // object.
  std::vector<double> operator()(double lc, double lt, double zt) const;

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

// We write using declarations so that we don't have to type the namespace name
// each time we use these names
using cosmosis::DataBlock;
using cubacpp::integration_results_v;

// We put the module label "modulelabel" in an anonymous namespace to make sure
// no other compilation unit can see it, and so that it has static lifetime.
// This means it is created when the library is loaded, and is destroyed when
// the program is shut down.
namespace {
  // We make the modulelabel const because it should never be modified.
  std::string const modulelabel("y1_analysis");
} // anonymous namespace

y1_analysis::y1_analysis(DataBlock& )
  : lc_lt()
{}

void
y1_analysis::set_sample(DataBlock& sample)
{
  // If we had a data member of type std::optional<X>, we would set the
  // value using std::optional::emplace(...) here. emplace takes a set
  // of arguments that it passes to the constructor of X.
  lc_lt.emplace(sample);
}

std::vector<double>
y1_analysis::operator()(double lc, double lt, double zt) const
{
  // For any data members of type std::optional<X>, we have to use operator*
  // to access the X object (as if we were dereferencing a pointer).
  std::vector<double> results(1);
  results[0] = (*lc_lt)(lc, lt, zt);
  return results;
}

//
void
y1_analysis::finalize_sample(cosmosis::DataBlock& sample,
                                  std::vector<integration_results_v> const& results)
  const
{
  // TODO: fix this to handle the whole vector of integration_results.
  std::vector<double> N_vals;
  std::vector<double> N_errors;
  std::vector<double> N_probs;
  std::vector<int> N_status;
  for (auto const &result : results){
      N_vals.push_back(result.value[0]);
      N_errors.push_back(result.error[0]);
      N_probs.push_back(result.prob[0]);
      N_status.push_back(result.status);
  } 
  sample.put_val(modulelabel, "N_vals", N_vals);
  sample.put_val(modulelabel, "N_errors", N_errors);
  sample.put_val(modulelabel, "N_probs", N_probs);
  sample.put_val(modulelabel, "N_status", N_status);
}

// The implementation of make_integration_volumes can be almost the same for
// any CosmoSISIntegrand-type class. Only the names and number of the parameters
// provided need to be changed. It is critical that the names be given in the
// order that correspond to the order of arguments in the class's function call
// operator. While the compiler can verify the number of arguments provided is
// correct, it can not verify that their order matches the order of arguments in
// the function call operator.
std::vector<y1_analysis::volume_t>
y1_analysis::make_integration_volumes(cosmosis::DataBlock& cfg)
{
  try{
  return y3_cluster::make_integration_volumes(cfg, modulelabel, "lc", "lt", "zt");
  } catch (std::exception const & ex){
	    std::cerr << ex.what();
             throw ;
      };
}


#include "CosmoSISIntegrationModule.hh"
#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "/cosmosis/cosmosis/datablock/entry.hh"
#include "/cosmosis/cosmosis/datablock/section.hh"
#include "/cosmosis/cosmosis/datablock/section_names.h"
#include <exception>
#include <iostream>

using y1_analysis_module =
  y3_cluster::CosmoSISIntegrationModule<y1_analysis>;

extern "C"
{
  void*
  setup(cosmosis::DataBlock* config)
  {
    try {
    return new y1_analysis_module(*config);
    } catch (std::exception const & ex){
	    std::cerr << ex.what();
             throw ;
      };
  }

  DATABLOCK_STATUS
  execute(cosmosis::DataBlock* sample, void* module)
  {
    try {
    auto mod = static_cast<y1_analysis_module*>(module);
    mod->execute(*sample);
    return DBS_SUCCESS;
    } catch (std::exception const & ex){
	    std::cerr << ex.what();
             throw ;
      };
  }

  int
  cleanup(void* module)
  {
    delete static_cast<y1_analysis_module*>(module);
    return 0;
  }
}


