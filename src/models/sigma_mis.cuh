#ifndef Y3_CLUSTER_SIGMA_MIS_CUH
#define Y3_CLUSTER_SIGMA_MIS_CUH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "common/cuda/Interp2D.cuh"
#include "models/ez.hh"
#include "utils/make_interp_2d.cuh"
#include "utils/primitives.hh"
#include "models/nfw_sigma_mis.cuh"

#include <algorithm>

namespace y3_cuda {
  class SIGMA_MIS {
  private:
    y3_cuda::NFW_SIGMA_MIS _nfw_sigma_mis;
    quad::Interp2D _sigma2;
    quad::Interp2D _bias;
    quad::Interp2D _sigma_crit_inv;
    double const _tau;

  public:
    size_t
    get_device_mem_footprint()
    {
      size_t size = 0;
      size += _nfw_sigma_mis.get_device_mem_footprint();
      size += _sigma2.get_device_mem_footprint();
      size += _bias.get_device_mem_footprint();
      size += _sigma_crit_inv.get_device_mem_footprint();
      return size;
    }

    SIGMA_MIS(quad::Interp2D const& nfw_sigma_mis,
              quad::Interp2D const& sigma2,
              quad::Interp2D const& bias,
              quad::Interp2D const& sigma_crit_inv)
      : _nfw_sigma_mis(nfw_sigma_mis), _sigma2(sigma2), _bias(bias), _sigma_crit_inv(sigma_crit_inv)
    {}

    using doubles = std::vector<double>;

    explicit SIGMA_MIS(cosmosis::DataBlock& sample)
      : _sigma2(make_Interp2D(sample,
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
      , _nfw_sigma_mis(y3_cuda::NFW_SIGMA_MIS())
      , _tau(sample.view<double>("cluster_abundance", "roffset_tau"))
    {}

    __device__ __host__ double
    operator()(double r, double lnM, double zt) const
    /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
    {
      // eval the mis-centered NFW
      // TODO: Implement the mass-concentration relation
      double const sig_mis = _nfw_sigma_mis.clamp(r, _tau, lnM);
      double const sig_2 = _bias.clamp(zt, lnM) * _sigma2.clamp(r, zt);
      double const sigc_inv = _sigma_crit_inv.clamp(zt, r);
      double const sig = std::max(sig_mis, sig_2);
      double const res = sig * sigc_inv;
      return res;
    }
  };
}

#endif
