#ifndef Y3_CLUSTER_KAPPA_MAX_CUH
#define Y3_CLUSTER_KAPPA_MAX_CUH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "cuda/pagani/quad/GPUquad/Interp2D.cuh"
#include "models/ez.hh"
#include "utils/make_interp_2d.cuh"
#include "utils/primitives.hh"

#include <algorithm>

namespace y3_cuda {
  class KAPPA_MAX {
  private:
    quad::Interp2D _sigma1;
    quad::Interp2D _sigma2;
    quad::Interp2D _bias;
    quad::Interp2D _sigma_crit_inv;

  public:
    size_t
    get_device_mem_footprint()
    {
      size_t size = 0;
      size += _sigma1.get_device_mem_footprint();
      size += _sigma2.get_device_mem_footprint();
      size += _bias.get_device_mem_footprint();
      return size;
    }

    KAPPA_MAX(quad::Interp2D const& sigma1,
              quad::Interp2D const& sigma2,
              quad::Interp2D const& bias,
              quad::Interp2D const& sigma_crit_inv)
      : _sigma1(sigma1), _sigma2(sigma2), _bias(bias), _sigma_crit_inv(sigma_crit_inv)
    {}

    using doubles = std::vector<double>;

    explicit KAPPA_MAX(cosmosis::DataBlock& sample)
      : _sigma1(make_Interp2D(sample,
                              "correlationFunction",
                              "r_sigma",
                              "lnM",
                              "Sigma_nfw"))
      , _sigma2(make_Interp2D(sample,
                              "correlationFunction",
                              "r_sigma",
                              "z",
                              "Sigma_hh"))
      , _bias(make_Interp2D(sample,
                            "correlationFunction",
                            "z",
                            "lnM",
                            "bias"))
      , _sigma_crit_inv(make_Interp2D(sample,
                              "correlationFunction",
                              "z",
                              "correlationFunction",
                              "r_sigma",
                              "sigmaCritInv",
                              "sigma_crit_inv"))                        
    {}

    __device__ __host__ double
    operator()(double r, double lnM, double zt) const
    /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
    {
      double const sig_1 = _sigma1.clamp(r, lnM);
      double const sig_2 = _bias.clamp(zt, lnM) * _sigma2.clamp(r, zt);
      double const sig = std::max(sig_1, sig_2);
      double const res = sig*_sigma_crit_inv(zt, r);
      return res;
    }
  };
}

#endif
