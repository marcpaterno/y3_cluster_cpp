#ifndef Y3_CLUSTER_DSIGMA_MISC_CUH
#define Y3_CLUSTER_DSIGMA_MISC_CUH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "utils/primitives.hh"
#include "models/ez.hh"
#include "models/nfw_dsigma_mis.cuh"

// ask Marc the new way to do this
#include "common/cuda/Interp2D.cuh"
#include "utils/cuda_interp_2d.cuh"
#include "utils/make_interp_2d.cuh"

#include <algorithm>
#include <fstream>
//  Mis-Centered Delta Sigma
//  Sigma_mis

namespace y3_cuda {
  class DSIGMA_MISC {
  private:
    y3_cuda::NFW_DSIGMA_MIS _nfw_dsigma_mis;
    quad::Interp2D _sigma2;
    quad::Interp2D _bias;
    quad::Interp2D _sigma_crit_inv;
    double const _tau;

  public:
    size_t
    get_device_mem_footprint()
    {
      size_t size = 0;
      // size += _nfw_dsigma_mis.get_device_mem_footprint();
      size += _sigma2.get_device_mem_footprint();
      size += _bias.get_device_mem_footprint();
      size += _sigma_crit_inv.get_device_mem_footprint();
      return size;
    }

    DSIGMA_MISC(y3_cuda::NFW_DSIGMA_MIS const& nfw_dsigma_mis,
                quad::Interp2D const& sigma2,
                quad::Interp2D const& bias,
                quad::Interp2D const& sigma_crit_inv,
                double const& tau
              )
      : _nfw_dsigma_mis(nfw_dsigma_mis), _sigma2(sigma2), 
        _bias(bias), _sigma_crit_inv(sigma_crit_inv), _tau(tau)
    {}

    using doubles = std::vector<double>;

    explicit DSIGMA_MISC(cosmosis::DataBlock& sample)
      : _nfw_dsigma_mis()
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
      // , _nfw_dsigma_mis(y3_cuda::NFW_DSIGMA_MIS(sample))
      , _tau(sample.view<double>("cluster_abundance", "roffset_tau"))
    {}

    __device__ __host__ double
    operator()(double r, double lnM, double zt) const
    /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
    {
      double const sig_2 = _bias.clamp(zt, lnM) * _sigma2.clamp(r, zt);
      double const sigc_inv = _sigma_crit_inv.clamp(zt, r);
      // eval the mis-centered NFW
      double const sig_mis_1 = _nfw_dsigma_mis(r, _tau, lnM);
      double const sig_mis = std::max(sig_mis_1, sig_2);
      double const res = sig_mis * sigc_inv;
      return res;
    }
  };
}

#endif
