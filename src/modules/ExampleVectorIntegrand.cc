#include "ExampleVectorIntegrand.hh"
#include "utils/make_integration_volumes.hh"

// We write using declarations so that we don't have to type the namespace name
// each time we use these names
using cosmosis::DataBlock;
using cubacpp::integration_results_v;

ExampleVectorIntegrand::ExampleVectorIntegrand(DataBlock& config)
  : radii_(config.view<std::vector<double>>(module_label(), "radii")), sigma_8_()
{}

void
ExampleVectorIntegrand::set_sample(DataBlock& sample)
{
  sigma_8_ = sample.view<double>("cosmological_parameters", "sigma_8");
  // If we had a data member of type std::optional<X>, we would set the
  // value using std::optional::emplace(...) here. emplace takes a set
  // of arguments that it passes to the constructor of X.
}

// This math is totally non-physical and stupid, but it uses all the values
// provided.
std::vector<double>
ExampleVectorIntegrand::operator()(double x, double y) const
{
  // For any data members of type std::optional<X>, we have to use operator*
  // to access the X object (as if we were dereferencing a pointer).
  std::vector<double> results(radii_.size());
  for (std::size_t i = 0; i != radii_.size(); ++i) {
    auto const delta = y - sigma_8_;
    results[i] = (x / radii_[i]) + (delta * delta);
  }
  return results;
}

//
void
ExampleVectorIntegrand::finalize_sample(cosmosis::DataBlock& sample,
                                  std::vector<integration_results_v> const& results)
  const
{
  // TODO: fix this to handle the whole vector of integration_results.
  sample.put_val(module_label(), "integral_vals", results[0].value);
  sample.put_val(module_label(), "integral_errors", results[0].error);
  sample.put_val(module_label(), "integral_probs", results[0].prob);
  sample.put_val(module_label(), "integral_status", results[0].status);
}


char const*
ExampleVectorIntegrand::module_label()
{
  return "example_vector_integrand";
}

// The implementation of make_integration_volumes can be almost the same for
// any CosmoSISIntegrand-type class. Only the names and number of the parameters
// provided need to be changed. It is critical that the names be given in the
// order that correspond to the order of arguments in the class's function call
// operator. While the compiler can verify the number of arguments provided is
// correct, it can not verify that their order matches the order of arguments in
// the function call operator.
std::vector<ExampleVectorIntegrand::volume_t>
ExampleVectorIntegrand::make_integration_volumes(cosmosis::DataBlock& cfg)
{
  return y3_cluster::make_integration_volumes_wall_of_numbers(cfg,
                                              ExampleVectorIntegrand::module_label(),
                                              "x",
                                              "y");
}
