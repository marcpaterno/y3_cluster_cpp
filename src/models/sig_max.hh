#ifndef Y3_CLUSTER_DEL_SIG_TOM_HH
#define Y3_CLUSTER_DEL_SIG_TOM_HH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "models/ez.hh"
#include "utils/interp_2d.hh"
#include "utils/make_interp_2d.hh"
#include "utils/primitives.hh"

#include <cmath>

namespace y3_cluster {
  class SIG_MAX {
  private:
    Interp2D _sigma1;
    Interp2D _sigma2;
    Interp2D _bias;

  public:
    SIG_MAX(Interp2D const& sigma1,
            Interp2D const& sigma2,
            Interp2D const& bias)
      : _sigma1(sigma1), _sigma2(sigma2), _bias(bias)
    {}

    using doubles = std::vector<double>;

    explicit SIG_MAX(cosmosis::DataBlock& sample)
      : _sigma1(make_Interp2D(sample, "deltasigma", "r_sigma_deltasigma", "lnM", "sigma_1"))
      , _sigma2(make_Interp2D(sample, "deltasigma", "r_sigma_deltasigma","matter_power_lin", "z","deltasigma","sigma_2"))
      , _bias(make_Interp2D(sample, "matter_power_lin", "z", "deltasigma", "lnM", "deltasigma", "bias"))
    {}

    double
    operator()(double r, double lnM, double zt) const
    /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
    {
      double _sig_1 = _sigma1.clamp(r, lnM);
      double _sig_2 = _bias.clamp(zt, lnM) * _sigma2.clamp(r, zt);
      double res = _sig_2;
      if (_sig_1 >= _sig_2) { res = _sig_1; }
      return res;
    }
  };
}

#endif
