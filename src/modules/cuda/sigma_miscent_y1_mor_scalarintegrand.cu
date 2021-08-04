#include "utils/cuda_module_macros.cuh"
#include "utils/make_cuda_integration_volumes.cuh"
#include "utils/make_grid_points.hh"

#include "cosmosis/datablock/datablock.hh"
#include "cubacpp/integration_result.hh"
#include "cudaPagani/quad/util/Volume.cuh"

#include "utils/make_interp_1d.cuh"

#include "models/dv_do_dz_t.cuh"
#include "models/hmf_t.cuh"
#include "models/int_lc_lt_des_t.cuh"
#include "models/int_zo_zt_des_t.cuh"
#include "models/lo_lc_t.cuh"
#include "models/mor_t.cuh"
#include "models/omega_z_des.cuh"
#include "models/roffset_t.cuh"
#include "models/sig_sum.cuh"
#
#include <iostream>
#include <optional>
#include <vector>

using namespace y3_cluster;
using cosmosis::DataBlock;
using cosmosis::ndarray;
using cubacpp::integration_result;

// This is a class that models the concept of
// "CosmoSISScalarIntegrand", and is thus suitable for use as the template
// parameter for the class template CosmoSISSICUDAModule.
//
// CAUTION: This is *very* preliminary. Don't hurt yourself.
//
class SigmaMiscentY1CUDAIntegrand {
public:
  using grid_t = y3_cluster::grid_t<3>;
  using grid_point_t = grid_t::value_type;

private:
  // We define the type alias volume_t to be the right dimensionality
  // of integration volume for our integrand. If we were to change the
  // number of arguments required by the function call operator (below),
  // we would need to also modify this type alias to keep consistent.
  using volume_t = quad::Volume<double, 6>;

  // State obtained from configuration. These things should be set in the
  // constructor.
  // <none in this example>

  // State obtained from each sample.
  // If there were a type X that did not have a default constructor,
  // we would use std::optional<X> as our data member.
  std::optional<INT_LC_LT_DES_t> lc_lt;
  std::optional<MOR_t> mor;
  std::optional<OMEGA_Z_DES> omega_z;
  std::optional<DV_DO_DZ_t> dv_do_dz;
  std::optional<HMF_t> hmf;
  std::optional<INT_ZO_ZT_DES_t> int_zo_zt;
  std::optional<ROFFSET_t> roffset;
  std::optional<LO_LC_t> lo_lc;
  std::optional<SIG_SUM> sigma;

  // State set for current 'bin' to be integrated.
  double zo_low_;
  double zo_high_;
  double radius_;

public:
  // Initialize my integrand object from the parameters read
  // from the relevant block in the CosmoSIS ini file.
  explicit SigmaMiscentY1CUDAIntegrand(cosmosis::DataBlock& config);

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
  __device__ double operator()(double lo,
                               double lc,
                               double zt,
                               double lnM,
                               double rmis,
                               double theta) const;

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

SigmaMiscentY1CUDAIntegrand::SigmaMiscentY1CUDAIntegrand(DataBlock&)
  : lc_lt()
  , mor()
  , omega_z()
  , dv_do_dz()
  , hmf()
  , int_zo_zt()
  , roffset()
  , lo_lc()
  , sigma()
  , zo_low_()
  , zo_high_()
  , radius_()
{}

void
SigmaMiscentY1CUDAIntegrand::set_sample(DataBlock& sample)
{
  // If we had a data member of type std::optional<X>, we would set the
  // value using std::optional::emplace(...) here. emplace takes a set
  // of arguments that it passes to the constructor of X.
  lc_lt.emplace(sample);
  mor.emplace(sample);
  dv_do_dz.emplace(sample);
  hmf.emplace(sample);
  omega_z.emplace(sample);
  roffset.emplace(sample);
  lo_lc.emplace(sample);
  sigma.emplace(sample);
}

void
SigmaMiscentY1CUDAIntegrand::set_grid_point(grid_point_t const& grid_point)
{
  radius_ = grid_point[2];
  zo_low_ = grid_point[0];
  zo_high_ = grid_point[1];
}

__device__  bool isbad(const char* name, double x) {
  if (isnan(x)) {
    printf("%s is nan\n", name);
    return true;
  }
  if (x == 0.0) {
    printf("%s is zero\n", name);
    return true;
  }
  if (abs(x) == std::numeric_limits<double>::infinity()) {
    printf("%s is infinite\n", name);
    return true;
  }
  return false;
}


__device__ double
SigmaMiscentY1CUDAIntegrand::operator()(double lo,
                                        double lc,
                                        double zt,
                                        double lnM,
                                        double rmis,
                                        double theta) const
{
  double const roffset_ = (*roffset)(rmis);
  if (isbad("roffset_", roffset_)) {
    printf("roffset bad: rmis=%g\n", rmis);
  }
  double const lo_lc_ = (*lo_lc)(lo, lc, rmis);
  if (isbad("lo_lc_", lo_lc_)) {
    printf("lo_lc bad: lo_lc_=%g  lo=%g lc=%g rmis=%g\n", lo_lc_, lo, lc, rmis);
  }
  double const mor_ = (*mor)(lc, lnM, zt);
  if (isbad("mor_", mor_)) {
    printf("mor bad: lo=%g lnM=%g zt=%g\n", lc, lnM, zt);
  }
  double const dv_do_dz_ = (*dv_do_dz)(zt);
  if (isbad("dv_do_dz_", dv_do_dz_)) {
    printf("dv_do_dz bad: zt=%g\n", zt);
  }
  double const hmf_ = (*hmf)(lnM, zt);
  if (isbad("hmf_", hmf_)) {
    printf("hmf bad: lnM=%g zt=%g\n", lnM, zt);
  }
  double const omega_z_ = (*omega_z)(zt);
  if (isbad("omega_z_", omega_z_)) {
    printf("omega_z bad: zt=%g\n", zt);
  }
  double common_term = roffset_ * lo_lc_ *
                       mor_ * dv_do_dz_ * hmf_ *
                       omega_z_ / 2.0 / 3.1415926535897;
  if (isbad("common_term", common_term)) {
    printf("common_term bad roffset_=%g lo_lc_=%g mor_=%g dv_do_dz_=%g hmf_=%g omega_z_=%g\n", roffset_, lo_lc_, mor_, dv_do_dz_, hmf_, omega_z_);
  }
  double scaled_Rmis = std::sqrt(radius_ * radius_ + rmis * rmis +
                                 2 * rmis * radius_ * std::cos(theta));
  double const sigma_ = (*sigma)(scaled_Rmis, lnM, zt);
  if (isbad("sigma_", sigma_)) {
    printf("sigma bad: radius_=%g rmis=%g theta=%g scaled_Rmis=%g lnM=%g zt=%g\n", radius_, rmis, theta, scaled_Rmis, lnM, zt);
  }
  double const int_zo_zt_ = (*int_zo_zt)(zo_low_, zo_high_, zt);
  if (isbad("int_zo_zt_", int_zo_zt_)) {
    printf("int_zo_zt bad: zo_low_=%g zo_high_=%g zt=%g\n", zo_low_, zo_high_, zt);
  }
  double const val = sigma_ * int_zo_zt_ * common_term;
  if (isbad("val", val)) {
      printf("operator(): generated a bad\n");
  }
  return val;
}

char const*
SigmaMiscentY1CUDAIntegrand::module_label()
{
  return "SigmaMiscentY1CUDAIntegrand";
}

std::vector<SigmaMiscentY1CUDAIntegrand::volume_t>
SigmaMiscentY1CUDAIntegrand::make_integration_volumes(cosmosis::DataBlock& cfg)
{
  return make_cuda_integration_volumes_wall_of_numbers(
    cfg,
    SigmaMiscentY1CUDAIntegrand::module_label(),
    "lo",
    "lc",
    "zt",
    "lnm",
    "rmis",
    "theta");
}

SigmaMiscentY1CUDAIntegrand::grid_t
SigmaMiscentY1CUDAIntegrand::make_grid_points(cosmosis::DataBlock& cfg)
{
  return y3_cluster::make_grid_points_cartesian_product(
    cfg,
    SigmaMiscentY1CUDAIntegrand::module_label(),
    "zo_low",
    "zo_high",
    "radii");
}

DEFINE_COSMOSIS_CUDA_INTEGRATION_MODULE(SigmaMiscentY1CUDAIntegrand)
