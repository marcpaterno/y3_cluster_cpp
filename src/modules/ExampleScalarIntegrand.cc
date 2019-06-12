#include "ExampleScalarIntegrand.hh"
#include "utils/make_integration_volumes.hh"

// We write using declarations so that we don't have to type the namespace name
// each time we use these names
using cosmosis::DataBlock;
using cubacpp::integration_result;

ExampleScalarIntegrand::ExampleScalarIntegrand(DataBlock& /* config */)
  : sigma_8_(), radius_()
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
ExampleScalarIntegrand::set_bin(bin_t const& bin)
{
  radius_ = bin[0];
}

// This math is totally non-physical and stupid, but it uses all the values
// provided.
double
ExampleScalarIntegrand::operator()(double x, double y) const
{
  // For any data members of type std::optional<X>, we have to use operator*
  // to access the X object (as if we were dereferencing a pointer).
  auto const delta = y - sigma_8_;
  return  (x / radius_) + (delta * delta);
}

//
void
ExampleScalarIntegrand::finalize_sample(cosmosis::DataBlock& sample,
                                  std::vector<integration_result> const& results)
  const
{
  // TODO: fix this to handle the whole vector of integration_results.
  sample.put_val(module_label(), "integral_vals", results[0].value);
  sample.put_val(module_label(), "integral_errors", results[0].error);
  sample.put_val(module_label(), "integral_probs", results[0].prob);
  sample.put_val(module_label(), "integral_status", results[0].status);
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
  return y3_cluster::make_integration_volumes(cfg,
                                              ExampleScalarIntegrand::module_label(),
                                              "x",
                                              "y");
}
