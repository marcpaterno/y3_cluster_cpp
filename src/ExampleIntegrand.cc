#include "ExampleIntegrand.hh"
#include "make_integration_volumes.hh"

// We write using declarations so that we don't have to type the namespace name
// each time we use DataBlock or make_integration_range.
using cosmosis::DataBlock;

// We put the module label "modulelabel" in an anonymous namespace to make sure
// no other compilation unit can see it, and so that it has static lifetime.
// This means it is created when the library is loaded, and is destroyed when
// the program is shut down.
namespace {
  // We make the modulelabel const because it should never be modified.
  std::string const modulelabel("example_integrand");
} // anonymous namespace

ExampleIntegrand::ExampleIntegrand(DataBlock& config)
  : radii_(config.view<std::vector<double>>(modulelabel, "radii")), sigma_8_()
{}

void
ExampleIntegrand::set_sample(DataBlock& sample)
{
  sigma_8_ = sample.view<double>("cosmological_parameters", "sigma_8");
  // If we had a data member of type std::optional<X>, we would set the
  // value using std::optional::emplace(...) here. emplace takes a set
  // of arguments that it passes to the constructor of X.
}

// This math is totally non-physical and stupid, but it uses all the values
// provided.
std::vector<double>
ExampleIntegrand::operator()(double x, double y) const
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

// The implementation of make_integration_volumes can be almost the same for
// any CosmoSISIntegrand-type class. Only the names and number of the parameters
// provided need to be changed. It is critical that the names be given in the
// order that correspond to the order of arguments in the class's function call
// operator. While the compiler can verify the number of arguments provided is
// correct, it can not verify that their order matches the order of arguments in
// the function call operator.
std::vector<ExampleIntegrand::volume_t>
ExampleIntegrand::make_integration_volumes(cosmosis::DataBlock& cfg)
{
  return y3_cluster_cpp::make_integration_volumes(cfg, modulelabel, "x", "y");
}
