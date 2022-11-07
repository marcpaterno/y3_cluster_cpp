#include "utils/cuda_module_macros.cuh"
#include "utils/datablock_reader.hh"
#include "utils/make_cuda_integration_volumes.cuh"
#include "utils/make_grid_points.hh"

#include "cubacpp/integration_result.hh"
#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/datablock_status.h"
#include "cuda/pagani/quad/util/Volume.cuh"

#include "models/dv_do_dz_t.cuh"
#include "models/hmf_t.cuh"
#include "models/int_zo_zt_des_t.cuh"
#include "models/lo_lc_t.cuh"
#include "models/mor_des_t.cuh"
#include "models/omega_z_des.cuh"
#include "models/roffset_t.cuh"

#include <iostream>
#include <optional>
#include <vector>

using cosmosis::DataBlock;
using cosmosis::ndarray;
using cubacpp::integration_result;

// summedNumbersCentY1GPU is a class that models the concept of
// "CosmoSISScalarIntegrand", and is thus suitable for use as the template
// parameter for the class template CosmoSISScalarIntegrationModule.
//
class summedNumbersCentY1GPU {
public:
  // Note we're slightly abusing the concept of the "grid" here. We are
  // really iterating over a single sequence of pairs, with the pairs
  // denoting the bottoms and tops of a set of zo bins.
  using grid_t = y3_cluster::grid_t<2>;
  using grid_point_t = grid_t::value_type;

private:
  // We define the type alias volume_t to be the right dimensionality
  // of integration volume for our integrand. If we were to change the
  // number of arguments required by the function call operator (below),
  // we would need to also modify this type alias to keep consistent.
  using volume_t = quad::Volume<double, 3>;

  // State obtained from configuration. These things should be set in the
  // constructor.
  // <none in this example>

  // State obtained from each sample.
  std::optional<y3_cuda::MOR_DES_t> mor;
  std::optional<y3_cuda::OMEGA_Z_DES> omega_z;
  std::optional<y3_cuda::DV_DO_DZ_t> dv_do_dz;
  std::optional<y3_cuda::HMF_t> hmf;
  std::optional<y3_cuda::INT_ZO_ZT_DES_t> int_zo_zt;

  // State set for current 'bin' to be integrated.
  double zo_low_;
  double zo_high_;


  // Should we use cartesian grid? If not, we use a wall-of-numbers.
  bool do_cartesian_product_of_bins_ = false;
  
  // In this integrand, the name "cartesian product" is a bit misleading; it
  // is used in case later development introduces another grid dimenstion over
  // which we'll operate.
  //
  // If we are using the cartesian grid, the relevant congifuration parameters
  // are zo_bottom (an array of arbitrary length)  and zo_bin_width (a array of
  // the same length as zo_bottom). These *two* parameters make a single axis
  // for the grid (the "zo bins").
  //
  // If we are using the wall-of-numbers configuration, the relevant
  // configuration parameters are zo_lo (an array of arbitrary length) and
  // zo_high. These *two* parameters make a single axis for the grid (the "zo
  // bins").

public:
  size_t
  get_device_mem_footprint()
  {
    size_t dev_size = 0;
    if (mor) dev_size += (*mor).get_device_mem_footprint();
    if (dv_do_dz) dev_size += (*dv_do_dz).get_device_mem_footprint();
    if (hmf) dev_size += (*hmf).get_device_mem_footprint();
    return dev_size;
  }

  // Initialize my integrand object from the parameters read
  // from the relevant block in the CosmoSIS ini file.
  explicit summedNumbersCentY1GPU(cosmosis::DataBlock& config);

  // Set any data members from values read from the current sample.
  // Do not attempt to copy the sample!.
  void set_sample(cosmosis::DataBlock& sample);

  // Set the data for the current bin.
  void set_grid_point(grid_point_t pt);

  // The function to be integrated. All arguments to this function must be of
  // type double, and there must be at least two of them (because our
  // integration routine does not work for functions of one variable). The
  // function is const because calling it does not change the state of the
  // object.
  __host__ __device__ double operator()(double lo, double zt, double lnM) const;

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

summedNumbersCentY1GPU::summedNumbersCentY1GPU(DataBlock& cfg)
{
  auto rc = 
    cfg.get_val(module_label(),
                "do_cartesian_product_of_bins",
                false,
                do_cartesian_product_of_bins_);
  if (rc != DBS_SUCCESS) {
    throw std::runtime_error("summedNumbersCentY1GPU failed to find do_cartesian_product_of_bins\n");
  }
}

void
summedNumbersCentY1GPU::set_sample(DataBlock& sample)
{
  // If we had a data member of type std::optional<X>, we would set the
  // value using std::optional::emplace(...) here. emplace takes a set
  // of arguments that it passes to the constructor of X.
  mor.emplace(sample);
  dv_do_dz.emplace(sample);
  hmf.emplace(sample);
  omega_z.emplace(sample);
}

void
summedNumbersCentY1GPU::set_grid_point(grid_point_t grid_point)
{
  zo_low_ = grid_point[0];
  zo_high_ = grid_point[1];
}

__host__ __device__ double
summedNumbersCentY1GPU::operator()(double lo,
                                            double zt,
                                            double lnM) const
{
  // For any data members of type std::optional<X>, we have to use operator*
  // to access the X object (as if we were dereferencing a pointer).
  double common_term =
    (*mor)(lo, lnM, zt) * (*dv_do_dz)(zt) * (*hmf)(lnM, zt) * (*omega_z)(zt);
  auto const val = (*int_zo_zt)(zo_low_, zo_high_, zt) * common_term;
  return val;
}

char const*
summedNumbersCentY1GPU::module_label()
{
  return "summedNumbersCentY1GPU";
}

// The implementation of make_integration_volumes can be almost the same for
// any CosmoSISIntegrand-type class. Only the names and number of the parameters
// provided need to be changed. It is critical that the names be given in the
// order that correspond to the order of arguments in the class's function call
// operator. While the compiler can verify the number of arguments provided is
// correct, it can not verify that their order matches the order of arguments in
// the function call operator.
std::vector<summedNumbersCentY1GPU::volume_t>
summedNumbersCentY1GPU::make_integration_volumes(
  cosmosis::DataBlock& cfg)
{
  return y3_cuda::make_integration_volumes_wall_of_numbers(
    cfg, summedNumbersCentY1GPU::module_label(), "lo", "zt", "lnm");
}

summedNumbersCentY1GPU::grid_t
summedNumbersCentY1GPU::make_grid_points(cosmosis::DataBlock& cfg)
{
  char const * const label =
    summedNumbersCentY1GPU::module_label();
  bool do_cartesian_product_of_bins = false;
  auto rc = cfg.get_val(label, "do_cartesian_product_of_bins", false, do_cartesian_product_of_bins);

  if (do_cartesian_product_of_bins)
  {
    auto const bottoms = get_vector_double(cfg, label, "zo_bottom");
    auto const widths= get_vector_double(cfg, label, "zo_bin_width");
    if (bottoms.size() != widths.size()) {
        throw std::runtime_error("zo_bottom and zo_bin_width must be the same length\n");
    }

    grid_t result;
    result.set_names(std::vector<std::string>{"zo_bottom", "zo_bin_width"});
    // This could probably be simplified by using the ranges library.
    for (std::size_t i = 0; i != bottoms.size(); ++i) {
      double bottom = bottoms[i];
      double top = bottoms[i] + widths[i];
      result.push_back({bottom, top});
    }
    return result;
  };

  return y3_cluster::make_grid_points_wall_of_numbers(
    cfg, summedNumbersCentY1GPU::module_label(), "zo_low", "zo_high");
}

DEFINE_COSMOSIS_CUDA_INTEGRATION_MODULE(summedNumbersCentY1GPU)
