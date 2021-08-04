#ifndef Y3_CLUSTER_CPP_HMF_T_CUH
#define Y3_CLUSTER_CPP_HMF_T_CUH

#include "cosmosis/datablock/datablock.hh"
#include "cudaPagani/quad/GPUquad/Interp2D.cuh"

#include <vector>

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
}

inline std::vector<double>
ndarray_z_to_vector(cosmosis::ndarray<double> const& a)
{
  return {a.begin(), a.end()};
}

class HMF_t {
public:
  HMF_t(quad::Interp2D const& nmz, double s, double q) : _nmz(nmz), _s(s), _q(q)
  {}

  using doubles = std::vector<double>;

  explicit HMF_t(cosmosis::DataBlock& sample)
    : _nmz(
        _adjust_to_log(sample, sample.view<doubles>("mass_function", "m_h")),
        sample.view<doubles>("mass_function", "z"),
        ndarray_z_to_vector(
          sample.view<cosmosis::ndarray<double>>("mass_function", "dndlnmh")))
    , _s(sample.view<double>("cluster_abundance", "hmf_s"))
    , _q(sample.view<double>("cluster_abundance", "hmf_q"))
  {}

  __device__ double
  operator()(double lnM, double zt) const
  {
    return _nmz.clamp(lnM, zt) *
           (_s * (lnM * 0.4342944819 - 13.8124426028) + _q);
    // 0.4342944819 is log(e)
  }

private:
  quad::Interp2D _nmz;
  double _s;
  double _q;
};

#endif
