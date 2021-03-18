#ifndef Y3_CLUSTER_DEL_SIG_T_HH
#define Y3_CLUSTER_DEL_SIG_T_HH

#include "cosmosis/datablock/datablock.hh"
#include "ez.hh"
#include "utils/datablock_reader.hh"
#include "utils/interp_2d.hh"
#include "utils/make_interp_2d.hh"
#include "utils/primitives.hh"

#include <cmath>
#include <memory>

namespace y3_cluster {
  class DEL_SIG_t {
  private:
    Interp2D _dsigma1;
    Interp2D _dsigma2;
    Interp2D _bias;

  public:
    DEL_SIG_t(Interp2D const& dsigma1,
              Interp2D const& dsigma2,
              Interp2D const& bias)
      : _dsigma1(dsigma1), _dsigma2(dsigma2), _bias(bias)
    {}

    using doubles = std::vector<double>;

    explicit DEL_SIG_t(cosmosis::DataBlock& sample)
      : _dsigma1(make_Interp2D(sample, "cluster_abundance", "del_sig_params_x1", "del_sig_params_y1", "del_sig_params_z1"))
      , _dsigma2(make_Interp2D(sample, "cluster_abundance", "del_sig_params_x2", "del_sig_params_y2", "del_sig_params_z2"))
      , _bias(make_Interp2D(sample, "cluster_abundance", "del_sig_params_x3", "del_sig_params_y3", "del_sig_params_z3"))
    {}

    double
    operator()(double r, double lnM, double zt) const
    /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
    {
      double del_sig_1 = _dsigma1(r, lnM);
      double del_sig_2 = _bias(zt, lnM) * _dsigma2(r, zt);
      if (del_sig_1 >= del_sig_2) {
        return (1. + zt) * (1. + zt) * (1. + zt) * del_sig_1;
      } else {
        return (1. + zt) * (1. + zt) * (1. + zt) * del_sig_2;
      }
    }
  };
}

#endif
