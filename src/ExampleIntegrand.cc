#include "ExampleIntegrand.hh"

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
  : radii_(config.view<double>(modulelabel, "radii"))
  , x_ranges_(make_integration_range(config, modulelabel, "x_range"))
  , y_ranges_(make_integration_range(config, modulelabel, "y_range"))
  , sigma_8_()
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
ExampleIntegrand::execute(double x, double y) const
{
  // For any data members of type std::optional<X>, we have to use operator*
  // to access the X object (as if we were dereferencing a pointer).
	std::vector<double> results(radii_.size());
  for (std::size_t i = 0; i != radii_.size(); ++i) {
	results[i] = (x/radii_[i]) + (y-sigma_8_);
  }
  return results;
}
