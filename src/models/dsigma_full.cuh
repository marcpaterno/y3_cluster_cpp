#ifndef Y3_CLUSTER_DSIGMA_FULL_CUH
#define Y3_CLUSTER_DSIGMA_FULL_CUH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "common/cuda/Interp2D.cuh"
#include "utils/make_interp_2d.cuh"
#include "utils/primitives.hh"
#include "models/ez.hh"
#include "models/nfw_dsigma_mis.cuh"

#include <algorithm>
#include <fstream>
// Joint Means Centered plus Mis-Centered 
// f_cen x Sigma + (1-f_cen) x Sigma_mis

namespace y3_cuda {
  class DSIGMA_FULL {
  private:
    y3_cuda::NFW_DSIGMA_MIS _nfw_dsigma_mis;
    quad::Interp2D _sigma1;
    quad::Interp2D _sigma2;
    quad::Interp2D _bias;
    quad::Interp2D _sigma_crit_inv;
    double const _fcen;
    double const _tau;

  public:
    size_t
    get_device_mem_footprint()
    {
      size_t size = 0;
      size += _sigma1.get_device_mem_footprint();
      // size += _nfw_dsigma_mis.get_device_mem_footprint();
      size += _sigma2.get_device_mem_footprint();
      size += _bias.get_device_mem_footprint();
      size += _sigma_crit_inv.get_device_mem_footprint();
      return size;
    }

    DSIGMA_FULL(
              quad::Interp2D const& sigma1,
              quad::Interp2D const& sigma2,
              quad::Interp2D const& bias,
              quad::Interp2D const& sigma_crit_inv,
              y3_cuda::NFW_DSIGMA_MIS const& nfw_dsigma_mis,
              double const& tau,
              double const& fcen
              )
      : _nfw_dsigma_mis(nfw_dsigma_mis), _sigma1(sigma1), _sigma2(sigma2), 
        _bias(bias), _sigma_crit_inv(sigma_crit_inv), _tau(tau), _fcen(fcen)
    {}

    using doubles = std::vector<double>;

    explicit DSIGMA_FULL(cosmosis::DataBlock& sample)
      : _sigma1(make_Interp2D(sample,
                              "correlationFunction",
                              "r_sigma",
                              "z",
                              "DSigma_nfw"))
      , _sigma2(make_Interp2D(sample,
                              "correlationFunction",
                              "r_sigma",
                              "z",
                              "DSigma_hh"))
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
      // , _nfw_dsigma_mis(y3_cuda::NFW_DSIGMA_MIS(sample))
      , _nfw_dsigma_mis()
      , _fcen(sample.view<double>("cluster_abundance", "fcen"))
      , _tau(sample.view<double>("cluster_abundance", "roffset_tau"))
    {}

    __device__ __host__ double
    operator()(double r, double lnM, double zt) const
    /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
    {
      double const sig_1 = _sigma1.clamp(r, lnM);
      // double const sig_1 = _nfw_dsigma_mis.clamp(r, 0.001, lnM);
      double const sig_2 = _bias.clamp(zt, lnM) * _sigma2.clamp(r, zt);
      double const sigc_inv = _sigma_crit_inv.clamp(zt, r);

      // eval the mis-centered NFW
      double const sig_mis_1 = _nfw_dsigma_mis(r, _tau, lnM);

      double const sig = std::max(sig_1, sig_2);
      double const sig_mis = std::max(sig_mis_1, sig_2);

      double const sig_joint = _fcen * sig + (1-_fcen) * sig_mis;
      double const res = sig_joint * sigc_inv;
      return res;
    }
  };
}

#endif
