#ifndef Y3_CLUSTER_SIGMA_MIS_JOINT_CUH
#define Y3_CLUSTER_SIGMA_MIS_JOINT_CUH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "common/cuda/Interp2D.cuh"
#include "utils/make_interp_2d.cuh"
#include "utils/primitives.hh"
#include "models/ez.hh"
#include "models/nfw_sigma_mis.cuh"

#include <algorithm>
// Joint Means Centered plus Mis-Centered 
// f_cen x Sigma + (1-f_cen) x Sigma_mis

namespace y3_cuda {
  class SIGMA_MIS_JOINT {
  private:
    y3_cuda::NFW_SIGMA_MIS _nfw_sigma_mis;
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
      size += _nfw_sigma_mis.get_device_mem_footprint();
      size += _sigma2.get_device_mem_footprint();
      size += _bias.get_device_mem_footprint();
      size += _sigma_crit_inv.get_device_mem_footprint();
      return size;
    }

    SIGMA_MIS_JOINT(quad::Interp2D const& sigma1,
              y3_cuda::NFW_SIGMA_MIS const& sigma2,
              quad::Interp2D const& sigma2,
              quad::Interp2D const& bias,
              quad::Interp2D const& sigma_crit_inv,
              double const& tau,
              double const& fcen)
      : _sigma1(sigma1), _nfw_sigma_mis(nfw_sigma_mis), _sigma2(sigma2), 
        _bias(bias), _sigma_crit_inv(sigma_crit_inv), _tau(tau), _fcen(fcen)
    {}

    using doubles = std::vector<double>;

    explicit SIGMA_MIS_JOINT(cosmosis::DataBlock& sample)
      : _sigma1(make_Interp2D(sample,
                              "correlationFunction",
                              "r_sigma",
                              "z",
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
      , _nfw_sigma_mis(y3_cuda::NFW_SIGMA_MIS(sample))
      , _fcen(sample.view<double>("cluster_abundance", "fcen"))
      , _tau(sample.view<double>("cluster_abundance" ,  "roffset_tau"))
    {}

    __device__ __host__ double
    operator()(double r, double lnM, double zt) const
    /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
    {
      double const sig_1 = _sigma1.clamp(r, lnM);
      // double const sig_1 = _nfw_sigma_mis.clamp(r, 0.001, lnM);
      double const sig_2 = _bias.clamp(zt, lnM) * _sigma2.clamp(r, zt);
      double const sigc_inv = _sigma_crit_inv.clamp(zt, r);

      // eval the mis-centered NFW
      double const sig_mis_1 = _nfw_sigma_mis.clamp(r, _tau, lnM);

      double const sig = std::max(sig_1, sig_2);
      double const sig_mis = std::max(sig_mis_1, sig_2);

      double const sig_joint = _fcen * sig + (1-_fcen) * sig_mis
      double const res = sig_joint * sigc_inv;
      return res;
    }
  };
}

#endif
