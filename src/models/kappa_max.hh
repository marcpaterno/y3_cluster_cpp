#ifndef Y3_CLUSTER_KAPPA_MAX_HH
#define Y3_CLUSTER_KAPPA_MAX_HH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "models/ez.hh"
#include "utils/interp_2d.hh"
#include "utils/make_interp_2d.hh"
#include "utils/primitives.hh"

#include <algorithm>

namespace y3_cluster {
  class KAPPA_MAX {
  private:
    Interp2D _sigma1;
    Interp2D _sigma2;
    Interp2D _bias;
    Interp2D _sigma_crit_inv;

  public:
    KAPPA_MAX(Interp2D const& sigma1,
              Interp2D const& sigma2,
              Interp2D const& bias,
              Interp2D const& sigma_crit_inv)
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
    {}

    double
    operator()(double r, double lnM, double zt) const
    /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
    {
      double const sig_1 = _sigma1.clamp(r, lnM);
      double const sig_2 = _bias.clamp(zt, lnM) * _sigma2.clamp(r, zt);
      double const sigc_inv = _sigma_crit_inv.clamp(zt, r);
      double const sig = std::max(sig_1, sig_2);
      double const res = sig * sigc_inv;
      return res;
    }
  };
}

#endif
