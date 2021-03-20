#ifndef Y3_CLUSTER_SIG_CRIT_INV
#define Y3_CLUSTER_SIG_CRIT_INV

#include "cosmosis/datablock/datablock.hh"
#include "utils/interp_1d.hh"
#include "utils/make_interp_1d.hh"
#include "utils/primitives.hh"

#include <vector>

namespace y3_cluster {

  class sigma_crit_inv {
  public:
    sigma_crit_inv(Interp1D const& da)
      : _da(da)
      // Constants are in primitives.hh
      // Units: Mpc/s (NO h)
      , _c(0.0)
      // Units: Mpc^3 / M_sol / s^2
      , _G(0.0)
    {}

    explicit sigma_crit_inv(cosmosis::DataBlock& sample)
      : sigma_crit_inv(make_Interp1D(sample, "distances", "z", "d_a"))
    {}

    double
    operator()(double zt, double zs) const
    {
      double _sig_crit_inv = 0;
      if (zs > zt) {
        // Units in Mpc, M_sol, s
        double const da_zt = _da(zt), // da_z needs to be in Mpc
          da_zs = _da(zs),            // da_z needs to be in Mpc
          da_zt_zs =
            da_zs - (1.0 + zt) / (1.0 + zs) * da_zt; // da_z needs to be in Mpc
        _sig_crit_inv = 4.0 * pi() * G() / c() / c() * da_zt * da_zt_zs / da_zs;
      }
      // Mpc^2 / M_sol
      return _sig_crit_inv;
    }

  private:
    Interp1D _da;
    double _c;
    double _G;
  };
}

#endif
