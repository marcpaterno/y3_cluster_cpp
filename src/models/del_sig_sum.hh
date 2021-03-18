#ifndef Y3_CLUSTER_DEL_SIG_TOM_HH
#define Y3_CLUSTER_DEL_SIG_TOM_HH

#include "ez.hh"
#include "utils/datablock_reader.hh"
#include "utils/interp_2d.hh"
#include "utils/primitives.hh"

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"

#include <cmath>

namespace y3_cluster {
  class DEL_SIG_TOM {
  private:
    Interp2D _dsigma1;
    Interp2D _dsigma2;
    Interp2D _bias;

  public:
    DEL_SIG_TOM(Interp2D const& dsigma1,
                Interp2D const& dsigma2,
                Interp2D const& bias)
      : _dsigma1(dsigma1), _dsigma2(dsigma2), _bias(bias)
    {}

    using doubles = std::vector<double>;

    explicit DEL_SIG_TOM(cosmosis::DataBlock& config)
      : _dsigma1(make_Interp2D(config, "deltasigma", "r_sigma_deltasigma", "lnM")),
      , _dsigma2(config.view<doubles>("deltasigma", "r_sigma_deltasigma"),
                 config.view<doubles>("matter_power_lin", "z"),
                 config.view<cosmosis::ndarray<double>>("deltasigma", "deltasigma_2"))
      , _bias(config.view<doubles>("matter_power_lin", "z"),
              config.view<doubles>("deltasigma", "lnM"),
              config.view<cosmosis::ndarray<double>>("deltasigma", "bias")))
    {}

    double
    operator()(double r, double lnM, double zt) const
    /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
    {
      double del_sig_1 = _dsigma1(r, lnM);
      double del_sig_2 = _bias(zt, lnM) * _dsigma2(r, zt);

      // NB: The 1e12 factor is to convert from the M_{sol} / kpc^2 of the input
      // to the M_{sol} / Mpc^2 we need
      // TODO: h factor?
      // TODO: Get deltasigma_2 working.
      // if (del_sig_1 >= del_sig_2) {
      return 1e12 * (1. + zt) * (1. + zt) * (1. + zt) * (del_sig_1 + del_sig_2);
      /*} else {
        return 1e12*(1.+zt)*(1.+zt)*(1.+zt)*del_sig_2;
      } */
    }
  };
}

#endif
