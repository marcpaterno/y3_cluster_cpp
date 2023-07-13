#ifndef Y3_CLUSTER_DSIGMA_PRJ_CUH
#define Y3_CLUSTER_DSIGMA_PRJ_CUH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "common/cuda/Interp2D.cuh"
#include "models/ez.hh"
#include "utils/make_interp_2d.cuh"
#include "utils/primitives.hh"
#include "models/nfw_dsigma_mis.cuh"

#include <algorithm>

namespace y3_cuda {
  class DSIGMA_PRJ {
  private:
    y3_cuda::NFW_DSIGMA_MIS _nfw_DSIGMA_MIS;
    quad::Interp2D _sigma2;
    quad::Interp2D _bias;
    quad::Interp2D _sigma_crit_inv;

  public:
    size_t
    get_device_mem_footprint()
    {
      size_t size = 0;
      size += _nfw_DSIGMA_MIS.get_device_mem_footprint();
      size += _sigma2.get_device_mem_footprint();
      size += _bias.get_device_mem_footprint();
      size += _sigma_crit_inv.get_device_mem_footprint();
      return size;
    }

    DSIGMA_PRJ(y3_cuda::NFW_DSIGMA_MIS const& nfw_DSIGMA_MIS,
              quad::Interp2D const& sigma2,
              quad::Interp2D const& bias,
              quad::Interp2D const& sigma_crit_inv)
      : _nfw_DSIGMA_MIS(nfw_DSIGMA_MIS), _sigma2(sigma2), _bias(bias), _sigma_crit_inv(sigma_crit_inv)
    {}

    using doubles = std::vector<double>;

    explicit DSIGMA_PRJ(cosmosis::DataBlock& sample)
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
      , _nfw_DSIGMA_MIS(y3_cuda::NFW_DSIGMA_MIS())
    {}

    __device__ __host__ double
    operator()(double r, double rmis, double lnM, double zt) const
    /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
    {
      // eval the mis-centered NFW; c=4 at z=0
      double const sig_mis = _nfw_DSIGMA_MIS.clamp(r, rmis, lnM);

      // piece-wise approximation
      // TODO: Eval the accuracy of the approx. for the 2h term
      double const r2 = r > rmis ? r : rmis;
      double const sig_2 = _bias.clamp(zt, lnM) * _sigma2.clamp(r2, zt);

      double const sig = std::max(sig_mis, sig_2);
      double const sigc_inv = _sigma_crit_inv.clamp(zt, r);
      double const res = sig * sigc_inv;
      return res;
    }
  };
}

#endif
