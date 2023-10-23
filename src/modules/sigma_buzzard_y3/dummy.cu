#include "utils/cuda_module_macros.cuh"
#include "utils/datablock_reader.hh"
#include "utils/make_cuda_integration_volumes.cuh"
#include "utils/make_grid_points.hh"

#include "cubacpp/integration_result.hh"
#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/datablock_status.h"
#include "common/cuda/Volume.cuh"

#include <vector>

using cosmosis::DataBlock;
using cosmosis::ndarray;
using cubacpp::integration_result;

// dummy is a class that is a very simple integrand and which uses
// none of the interpolation routines in y3_cluster_cpp. It is intended
// to help identify the source of the CPU memory exhaustion we are
// encountering.

class dummy {
public:
  using grid_t = y3_cluster::grid_t<1>;
  using grid_point_t = grid_t::value_type;

private:
  // We define the type alias volume_t to be the right dimensionality
  // of integration volume for our integrand. If we were to change the
  // number of arguments required by the function call operator (below),
  // we would need to also modify this type alias to keep consistent.
  using volume_t = quad::Volume<double, 2>;

    // State set for current 'bin' to be integrated.
  double a_;

  // Should we use cartesian grid? If not, we use a wall-of-numbers.
  bool do_cartesian_product_of_bins_ = false;
  
  public:
  size_t
  get_device_mem_footprint()
  {
    // Not sure what to return here. This is a best guess.
    return 8; // the size of our one double.
  }

  explicit dummy(cosmosis::DataBlock& cfg);
  void set_sample(cosmosis::DataBlock& sample);

  // Set the data for the current bin.
  void set_grid_point(grid_point_t pt);

  // The function to be integrated. All arguments to this function must be of
  // type double, and there must be at least two of them (because our
  // integration routine does not work for functions of one variable). The
  // function is const because calling it does not change the state of the
  // object.
  __host__ __device__ double operator()(double x, double y) const;

  // module_label() is a non-member (static) function that returns the label for
  // this module. The name this returns
  // is the name that must be used in the 'ini file' for configuring the module
  // made with this class.
  // We return char const* rather than std::string to avoid some needless memory
  // allocations.
  static char const* module_label();

  // The following non-member (static) function creates a vector of integration
  // volumes (the type alias defined above) based on the parameters read from
  // the configuration block for the module.
  static std::vector<volume_t> make_integration_volumes(
          cosmosis::DataBlock& cfg);

  // The following non-member (static) function creates a vector of grid points
  // on which the integration results are to be evaluated, based on parameters
  // read from the configuration block for the module.
  static grid_t make_grid_points(cosmosis::DataBlock& cfg);
};

// We write using declarations so that we don't have to type the namespace name
// each time we use these names
using cosmosis::DataBlock;
using cubacpp::integration_result;

dummy::dummy(DataBlock& cfg)
{
  auto rc = 
    cfg.get_val(module_label(),
                "do_cartesian_product_of_bins",
                false,
                do_cartesian_product_of_bins_);
  if (rc != DBS_SUCCESS) {
    throw std::runtime_error("dummy failed to find do_cartesian_product_of_bins\n");
  }
}

void
dummy::set_sample(DataBlock&)
{ }

void
dummy::set_grid_point(grid_point_t grid_point)
{
  a_ = grid_point[0];
}

__host__ __device__ double
dummy::operator()(double x, double y) const
{
  double const val = a_ * x + y*y;
  return val;
}

// string must match section block in pipeline.ini file
char const*
dummy::module_label()
{
  return "dummy";
}

// The implementation of make_integration_volumes can be almost the same for
// any CosmoSISIntegrand-type class. Only the names and number of the parameters
// provided need to be changed. It is critical that the names be given in the
// order that correspond to the order of arguments in the class's function call
// operator. While the compiler can verify the number of arguments provided is
// correct, it can not verify that their order matches the order of arguments in
// the function call operator.
std::vector<dummy::volume_t>
dummy::make_integration_volumes(
  cosmosis::DataBlock& cfg)
{
  return y3_cuda::make_integration_volumes_wall_of_numbers(
    cfg, 
    dummy::module_label(), 
    "x",
    "y");
}

dummy::grid_t
dummy::make_grid_points(cosmosis::DataBlock& cfg)
{
  char const * const label =
    dummy::module_label();

  return y3_cluster::make_grid_points_wall_of_numbers(
    cfg, 
    dummy::module_label(), 
    "a");
}

DEFINE_COSMOSIS_CUDA_INTEGRATION_MODULE(dummy)
