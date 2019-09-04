#include "ExampleScalarIntegrand.hh"
#include "utils/make_integration_volumes.hh"
#include "utils/make_grid_points.hh"
#include "/cosmosis/cosmosis/datablock/ndarray.hh"

// We write using declarations so that we don't have to type the namespace name
// each time we use these names
using cosmosis::DataBlock;
using cosmosis::ndarray;
using cubacpp::integration_result;

ExampleScalarIntegrand::ExampleScalarIntegrand(DataBlock&)
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
ExampleScalarIntegrand::set_grid_point(grid_point_t const& grid_point)
{
  radius_ = grid_point[0];
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
                                        std::vector<grid_point_t> const& grid_points,
                                        std::size_t nvolumes,
                                        std::vector<integration_result> const& res)
  const
{
  auto ngrid_points = grid_points.size();
  auto nresults = nvolumes * ngrid_points;

  // Create ndarray to give a view into the 'res' vector.
  ndarray<integration_result> results(res, {nvolumes, ngrid_points});

  // Create storage buffers for ndarrays to be inserted into the sample, and then
  // create the ndarrays.
  std::vector<double> vals_buffer(nresults);
  ndarray<double> vals(vals_buffer, {nvolumes, ngrid_points});
  
  std::vector<double> errors_buffer(nresults);
  ndarray<double> errors(errors_buffer, {nvolumes, ngrid_points});

  std::vector<double> probs_buffer(nresults);
  ndarray<double> probs(probs_buffer, {nvolumes, ngrid_points});

  std::vector<int> statuses_buffer(nresults);
  ndarray<int> statuses(statuses_buffer, {nvolumes, ngrid_points});

  for (std::size_t ivol = 0; ivol != nvolumes; ++ivol) {
    for (std::size_t jgp = 0; jgp != ngrid_points; ++jgp) {
      vals(ivol, jgp) = results(ivol, jgp).value;
      errors(ivol, jgp) = results(ivol, jgp).error;
      probs(ivol, jgp) = results(ivol, jgp).prob;
      statuses(ivol, jgp) = results(ivol, jgp).status;
    }
  }

  sample.put_val(module_label(), "integral_vals", vals);
  sample.put_val(module_label(), "integral_errors", errors);
  sample.put_val(module_label(), "integral_probs", probs);
  sample.put_val(module_label(), "integral_status", statuses);
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

std::vector<ExampleScalarIntegrand::grid_point_t>
ExampleScalarIntegrand::make_grid_points(cosmosis::DataBlock& cfg)
{
  return y3_cluster::make_grid_points(cfg,
                                      ExampleScalarIntegrand::module_label(),
                                      "radii");
}

