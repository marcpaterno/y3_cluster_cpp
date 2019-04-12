#include "ExampleIntegrand.hh"

// We write using declarations so that we don't have to type the namespace name
// each time we use DataBlock or make_integration_range.
using cosmosis::DataBlock;
using y3_cluster::make_integration_range;

// We put the module label "modulelabel" in an anonymous namespace to make sure
// no other compilation unit can see it, and so that it has static lifetime.
// This means it is created when the library is loaded, and is destroyed when
// the program is shut down.
namespace {
  // We make the modulelabel const because it should never be modified.
  std::string const modulelabel("example_integrand");
} // namespace

ExampleIntegrand::ExampleIntegrand(DataBlock& config)
  : radii_(config.view<double>(modulelabel, "radii"))
  , lnM_ranges_(make_integration_range(config, modulelabel, "lnM_range"))
  , z_ranges_(make_integration_range(config, modulelabel, "z_range"))
  , sigma_8_()
{}

void
ExampleIntegrand::set_sample(DataBlock& sample)
{
  sigma_8_ = sample.view<double>("cosmological_parameters", "sigma_8");
}

// This math is totally non-physical and stupid, but it uses all the values
// provided.
std::vector<double>
ExampleIntegrand::execute(double lnM, double zt) const
{
  // Note that because the type of hmf_ is std::optional<HMF_t>, not HMF_t, we
  // must use operator* on it to obtain the callable thing.
  return (*hmf_)(lnM + R_, zt + sigma_8_);
}
