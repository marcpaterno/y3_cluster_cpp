#ifndef Y3_CLUSTER_DSIGMA_MISC_HH
#define Y3_CLUSTER_DSIGMA_MISC_HH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "utils/primitives.hh"
#include "models/ez.hh"
#include "models/nfw_dsigma_mis.hh"

#include "utils/interp_2d.hh"
#include "utils/make_interp_2d.hh"

#include <algorithm>
#include <fstream>
//  Mis-Centered Delta Sigma
//  Sigma_mis

namespace y3_cluster {
  class DSIGMA_MISC {
  private:
    y3_cluster::NFW_DSIGMA_MIS _nfw_dsigma_mis;
    Interp2D _sigma2;
    Interp2D _bias;
    Interp2D _sigma_crit_inv;
    double const _tau;

  public:
    DSIGMA_MISC(y3_cluster::NFW_DSIGMA_MIS const& nfw_dsigma_mis,
                Interp2D const& sigma2,
                Interp2D const& bias,
                Interp2D const& sigma_crit_inv,
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
                              "DSigma_hh"))
      , _bias(make_Interp2D(sample,
                            "correlationFunction",
                            "lnM",
                            "z",
                            "bias"))
      , _sigma_crit_inv(make_Interp2D(sample,
                              "correlationFunction",
                              "r_sigma",
                              "correlationFunction",
                              "z",
                              "sigmaCritInv",
                              "sigma_crit_inv"))
      // , _nfw_dsigma_mis(y3_cuda::NFW_DSIGMA_MIS(sample))
      , _tau(sample.view<double>("cluster_abundance", "roffset_tau"))
    {}

    double
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
