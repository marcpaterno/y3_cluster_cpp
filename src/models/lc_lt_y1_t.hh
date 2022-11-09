#ifndef Y3_CLUSTER_LC_LT_Y1_T_HH
#define Y3_CLUSTER_LC_LT_Y1_T_HH

#include "cosmosis/datablock/datablock.hh"
#include "utils/interp_2d.hh"
#include "utils/primitives.hh"

// Implementing Costanzi, Rozo, Rykoff et al 2018 equation 6
// which is a P(\lambda^obs | \lambda^true, z) form

#include <cmath>

namespace y3_cluster {
  struct LC_LT_Y1_t {

    Interp2D tau_interp;
    Interp2D mu_interp;
    Interp2D sigma_interp;
    Interp2D fmsk_interp;
    Interp2D fprj_interp;

    LC_LT_Y1_t();
    // LC_LT_Y1_t(const cosmosis::DataBlock&) {TODO}

    double
    operator()(double lamcent, double lamtrue, double ztrue) const
    {
      // Values of each fit parameter at lamtrue, ztrue
      const auto tau = tau_interp(lamtrue, ztrue);
      const auto mu = mu_interp(lamtrue, ztrue);
      const auto sigma = sigma_interp(lamtrue, ztrue);
      const auto fmsk = fmsk_interp(lamtrue, ztrue);
      const auto fprj = fprj_interp(lamtrue, ztrue);

      // Some repeated expressions
      const auto exptau =
        std::exp(tau * (2.0 * mu + tau * sigma * sigma - 2.0 * lamcent) / 2.0);
      const auto root_two_sigma = std::sqrt(2.0) * sigma;
      const auto mu_tau_sig_sqr = mu + tau * sigma * sigma;

      // Helper function for common pattern
      const auto erfc_scaled = [root_two_sigma](double a, double b) {
        return std::erfc((a - b) / root_two_sigma);
      };

      return (1.0 - fmsk) * (1.0 - fprj) *
               y3_cluster::gaussian(lamcent, mu, sigma) +
             0.5 * ((1.0 - fmsk) * fprj * tau + fmsk * fprj / lamtrue) *
               exptau * erfc_scaled(mu_tau_sig_sqr, lamcent) +
             0.5 * fmsk / lamtrue *
               (erfc_scaled(lamcent, mu) - erfc_scaled(lamcent + lamtrue, mu)) -
             0.5 * fmsk * fprj / lamtrue *
               (std::exp(-tau * lamtrue) * exptau *
                erfc_scaled(mu_tau_sig_sqr, lamcent + lamtrue));
    }
  };
}

#endif
