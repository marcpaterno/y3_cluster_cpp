#include "utils/cuda_module_macros.cuh"
#include "utils/make_cuda_integration_volumes.cuh"
#include "utils/make_grid_points.hh"

#include "cosmosis/datablock/datablock.hh"
#include "cubacpp/integration_result.hh"
#include "cudaPagani/quad/util/Volume.cuh"

#include "models/dv_do_dz_t.cuh"
#include "models/hmf_t.cuh"
#include "models/int_lc_lt_des_t.cuh"
#include "models/int_zo_zt_des_t.cuh"
#include "models/mor_t.cuh"
#include "models/omega_z_des.cuh"
#include "models/roffset_t.cuh"
#include "models/sig_max.cuh"

#include <iostream>
#include <optional>
#include <vector>

using cosmosis::DataBlock;
using cosmosis::ndarray;
using cubacpp::integration_result;

// This is a class that models the concept of
// "CosmoSISCUDAScalarIntegrand", and is thus suitable for use as the template
// parameter for the class template CosmoSISScalarIntegrationModule.
//
class SigmaCentY1CUDAIntegrand {
public:
  using grid_t = y3_cluster::grid_t<3>;
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
  std::optional<y3_cuda::INT_LC_LT_DES_t> lc_lt;
  std::optional<y3_cuda::MOR_t> mor;
  std::optional<y3_cuda::OMEGA_Z_DES> omega_z;
  std::optional<y3_cuda::DV_DO_DZ_t> dv_do_dz;
  std::optional<y3_cuda::HMF_t> hmf;
  std::optional<y3_cuda::INT_ZO_ZT_DES_t> int_zo_zt;
  std::optional<y3_cuda::ROFFSET_t> roffset;
  std::optional<y3_cuda::SIG_MAX> sigma;

  // State set for current 'bin' to be integrated.
  double zo_low_;
  double zo_high_;
  double radius_;

public:
  // Initialize my integrand object from the parameters read
  // from the relevant block in the CosmoSIS ini file.
  explicit SigmaCentY1CUDAIntegrand(cosmosis::DataBlock& config);

  // Set any data members from values read from the current sample.
  // Do not attempt to copy the sample!.
  void set_sample(cosmosis::DataBlock& sample);

  // Set the data for the current bin.
  void set_grid_point(grid_point_t const& pt);

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

SigmaCentY1CUDAIntegrand::SigmaCentY1CUDAIntegrand(DataBlock&)
  : lc_lt()
  , mor()
  , omega_z()
  , dv_do_dz()
  , hmf()
  , int_zo_zt()
  , sigma()
  , zo_low_()
  , zo_high_()
  , radius_()
{}

void
SigmaCentY1CUDAIntegrand::set_sample(DataBlock& sample)
{
  // If we had a data member of type std::optional<X>, we would set the
  // value using std::optional::emplace(...) here. emplace takes a set
  // of arguments that it passes to the constructor of X.
  lc_lt.emplace(sample);
  mor.emplace(sample);
  dv_do_dz.emplace(sample);
  hmf.emplace(sample);
  omega_z.emplace(sample);
  sigma.emplace(sample);
}

void
SigmaCentY1CUDAIntegrand::set_grid_point(grid_point_t const& grid_point)
{
  radius_ = grid_point[2];
  zo_low_ = grid_point[0];
  zo_high_ = grid_point[1];
}

__host__ __device__ double
SigmaCentY1CUDAIntegrand::operator()(double lo,
                                           double zt,
                                           double lnM) const
{
  // For any data members of type std::optional<X>, we have to use operator*
  // to access the X object (as if we were dereferencing a pointer).
  double common_term =
    (*mor)(lo, lnM, zt) * (*dv_do_dz)(zt) * (*hmf)(lnM, zt) * (*omega_z)(zt);
  auto const val = (*sigma)(radius_, lnM, zt) *
                   (*int_zo_zt)(zo_low_, zo_high_, zt) * common_term;
  return val;
}

char const*
SigmaCentY1CUDAIntegrand::module_label()
{
  return "SigmaCentY1CUDAIntegrand";
}

// The implementation of make_integration_volumes can be almost the same for
// any CosmoSISIntegrand-type class. Only the names and number of the parameters
// provided need to be changed. It is critical that the names be given in the
// order that correspond to the order of arguments in the class's function call
// operator. While the compiler can verify the number of arguments provided is
// correct, it can not verify that their order matches the order of arguments in
// the function call operator.
std::vector<SigmaCentY1CUDAIntegrand::volume_t>
SigmaCentY1CUDAIntegrand::make_integration_volumes(
  cosmosis::DataBlock& cfg)
{
  return y3_cuda::make_integration_volumes_wall_of_numbers(
    cfg, SigmaCentY1CUDAIntegrand::module_label(), "lo", "zt", "lnm");
}

SigmaCentY1CUDAIntegrand::grid_t
SigmaCentY1CUDAIntegrand::make_grid_points(cosmosis::DataBlock& cfg)
{
  return y3_cluster::make_grid_points_wall_of_numbers(
    cfg,
    SigmaCentY1CUDAIntegrand::module_label(),
    "zo_low",
    "zo_high",
    "radii");
}

DEFINE_COSMOSIS_CUDA_INTEGRATION_MODULE(SigmaCentY1CUDAIntegrand)
