#ifndef Y3_CLUSTER_DEL_SIG_TOM_HH
#define Y3_CLUSTER_DEL_SIG_TOM_HH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "models/ez.hh"
#include "utils/datablock_reader.hh"
#include "utils/interp_2d.hh"
#include "utils/primitives.hh"

#include <cmath>
#include <memory>

namespace y3_cluster {
  class XI_MAX {
  private:
    std::shared_ptr<Interp2D const> _sigma1;
    std::shared_ptr<Interp2D const> _sigma2;
    std::shared_ptr<Interp2D const> _bias;
    double _om;

  public:
    XI_MAX(std::shared_ptr<Interp2D const> sigma1,
           std::shared_ptr<Interp2D const> sigma2,
           std::shared_ptr<Interp2D const> bias,
           double om)
      : _sigma1(sigma1), _sigma2(sigma2), _bias(bias), _om(om)
    {}

    using doubles = std::vector<double>;

    explicit XI_MAX(cosmosis::DataBlock& sample)
      : _sigma1(std::make_shared<Interp2D const>(
          get_datablock<doubles>(sample, "deltasigma", "r_xi"),
          get_datablock<doubles>(sample, "deltasigma", "lnM"),
          get_datablock<cosmosis::ndarray<double>>(sample,
                                                   "deltasigma",
                                                   "xi_1")))
      , _sigma2(std::make_shared<Interp2D const>(
          get_datablock<doubles>(sample, "deltasigma", "r_xi"),
          get_datablock<doubles>(sample, "matter_power_lin", "z"),
          get_datablock<cosmosis::ndarray<double>>(sample,
                                                   "deltasigma",
                                                   "xi_2")))
      , _bias(std::make_shared<Interp2D const>(
          get_datablock<doubles>(sample, "matter_power_lin", "z"),
          get_datablock<doubles>(sample, "deltasigma", "lnM"),
          get_datablock<cosmosis::ndarray<double>>(sample,
                                                   "deltasigma",
                                                   "bias")))
      , _om(get_datablock<double>(sample, "cosmological_parameters", "omega_m"))
    {}

    double
    operator()(double r, double lnM, double zt) const
    /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
    {
      double _sig_1 = _sigma1->eval(r, lnM);
      double _sig_2 = _bias->eval(zt, lnM) * _sigma2->eval(r, zt);
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
