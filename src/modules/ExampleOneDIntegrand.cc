#include "ExampleOneDIntegrand.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "utils/make_grid_points.hh"

// Constructor. Anything that need only be done once should be done here.
ExampleOneDIntegrand::ExampleOneDIntegrand(cosmosis::DataBlock& cfg)
  : b_(cfg.view<double>(module_label(), "b"))
  , c_(), a_()
{}

// The name of the module to use for configuration options and to
// save outputs into the datablock
char const*
ExampleOneDIntegrand::module_label()
{
  return "example_oned_integrand";
}

// The name of the integration variable to reference in the config file
char const*
ExampleOneDIntegrand::integration_variable()
{
  return "x";
}

// This is our integrand. In this example, we use f(x) = a*x + b + c
// For any data members of type std::optional<X>, we have to use operator*
// to access the X object (as if we were dereferencing a pointer).
double
ExampleOneDIntegrand::operator()(double x) const
{
  return a_*x + b_ + c_;
}

// Update anything that depends on the MCMC sample
// If we had a data member of type std::optional<X>, we would set the
// value using std::optional::emplace(...) here. emplace takes a set
// of arguments that it passes to the constructor of X.
void
ExampleOneDIntegrand::set_sample(cosmosis::DataBlock& sample)
{
  c_ = sample.view<double>(module_label(), "c");
}

// Update the grid points
void
ExampleOneDIntegrand::set_grid_point(grid_point_t const& grid_point)
{
  a_ = grid_point[0];
}

// Load the grid points from the config file. It is important that the
// order of the strings here is the same order that is used in the
// function set_grid_point
std::vector<ExampleOneDIntegrand::grid_point_t>
ExampleOneDIntegrand::make_grid_points(cosmosis::DataBlock& cfg)
{
  return y3_cluster::make_grid_points_wall_of_numbers(
    cfg, module_label(), "a");
}

// --------------------------------------------------------------------------
// ------- Function definitions below will be unchanged for most ------------
// --------------------------------------------------------------------------
// Load the integration ranges. This should not need to be adjusted.
std::vector<std::array<double, 2>>
ExampleOneDIntegrand::make_integration_volumes(cosmosis::DataBlock& cfg)
{
  auto modulelabel = module_label();
  std::string intvarname(integration_variable());
  std::string lowname = intvarname + "_low";
  std::string highname = intvarname + "_high";
  auto lows = get_vector_double(cfg, modulelabel, lowname);
  auto highs = get_vector_double(cfg, modulelabel, highname);

  std::vector<std::array<double, 2>> volumes;
  for (std::size_t i = 0; i != lows.size(); ++i)
    volumes.push_back(std::array<double, 2> { lows[i], highs[i] });

  return volumes;
}

