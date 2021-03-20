#ifndef Y3_CLUSTER_DEL_SIG_TOM_HH
#define Y3_CLUSTER_DEL_SIG_TOM_HH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "models/ez.hh"
#include "utils/interp_2d.hh"
#include "utils/make_interp_2d.hh"
#include "utils/primitives.hh"

#include <cmath>
#include <memory>

namespace y3_cluster {
  class XI_MAX {
  private:
    Interp2D _sigma1;
    Interp2D _sigma2;
    Interp2D _bias;
    double _om;

  public:
    XI_MAX(Interp2D const& sigma1,
           Interp2D const& sigma2,
           Interp2D const& bias,
           double om)
      : _sigma1(sigma1), _sigma2(sigma2), _bias(bias), _om(om)
    {}

    using doubles = std::vector<double>;

    explicit XI_MAX(cosmosis::DataBlock& sample)
      : _sigma1(make_Interp2D(sample, "deltasigma", "r_xi", "lnM", "xi_1"))
      , _sigma2(make_Interp2D(sample,
                              "deltasigma",
                              "r_xi",
                              "matter_power_lin",
                              "z",
                              "deltasigma",
                              "xi_2"))
      , _bias(make_Interp2D(sample,
                            "matter_power_lin",
                            "z",
                            "deltasigma",
                            "lnM",
                            "deltasigma",
                            "bias"))
      , _om(sample.view<double>("cosmological_parameters", "omega_m"))
    {}

    double
    operator()(double r, double lnM, double zt) const
    /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
    {
      double _sig_1 = _sigma1(r, lnM);
      double _sig_2 = _bias(zt, lnM) * _sigma2(r, zt);
      // TODO: h factor?
      double res = _sig_2;
      if (_sig_1 >= _sig_2) { res = _sig_1; }
      // return (1.+zt)*(1.+zt)*(1.+zt)*(_sig_1+_sig_2);
      return (res)*2.77536627E11 * _om / 1.0E10;
      /*} else {
        return 1e12*(1.+zt)*(1.+zt)*(1.+zt)*del_sig_2;
      } */
    }
  };
}

#endif
