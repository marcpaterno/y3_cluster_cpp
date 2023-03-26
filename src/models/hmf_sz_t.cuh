#ifndef Y3_CLUSTER_CPP_HMF_T_CUH
#define Y3_CLUSTER_CPP_HMF_T_CUH

#include "cosmosis/datablock/datablock.hh"
#include "common/cuda/Interp2D.cuh"

#include <vector>

// Use tinker halo-mass-function without s,q parametrization
// 

namespace {
  inline std::vector<double>
  _adjust_to_log(cosmosis::DataBlock& db, const std::vector<double>& masses)
  {
    std::vector<double> output = masses;
    double const omega_m =
      db.view<double>("cosmological_parameters", "omega_M");
    double const omega_mu =
      db.view<double>("cosmological_parameters", "omega_nu");
    for (auto& x : output) x = log(x * (omega_m - omega_mu));
    return output;
  }

  inline std::vector<double>
  ndarray_z_to_vector(cosmosis::ndarray<double> const& a)
  {
    return {a.begin(), a.end()};
  }
}

namespace y3_cuda {
  class HMF_t {
  public:
    HMF_t(quad::Interp2D const& nmz)
      : _nmz(nmz)
    {}

    using doubles = std::vector<double>;

    size_t
    get_device_mem_footprint()
    {
      return _nmz.get_device_mem_footprint();
    }

    explicit HMF_t(cosmosis::DataBlock& sample)
      : _nmz(
          _adjust_to_log(sample, sample.view<doubles>("mass_function", "m_h")),
          sample.view<doubles>("mass_function", "z"),
          ndarray_z_to_vector(
            sample.view<cosmosis::ndarray<double>>("mass_function", "dndlnmh")))
    {}

    __device__ __host__ double
    operator()(double lnM, double zt) const
    {
      return _nmz.clamp(lnM, zt)
    }

  private:
    quad::Interp2D _nmz;
  };
}

#endif
