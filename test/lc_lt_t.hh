#ifndef Y3_CLUSTER_LC_LT_T_HH
#define Y3_CLUSTER_LC_LT_T_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "test/interp_2d.hh"
#include "test/primitives.hh"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

namespace y3_cluster {
  struct LC_LT_t {

    static Interp2D const tau_interp;
    static Interp2D const mu_interp;
    static Interp2D const sigma_interp;
    static Interp2D const fmsk_interp;
    static Interp2D const fprj_interp;

    double
    operator()(double lc, double lt, double zt) const
    {
      const auto tau = tau_interp(lt, zt);
      const auto mu = mu_interp(lt, zt);
      const auto sigma = sigma_interp(lt, zt);
      const auto fmsk = fmsk_interp(lt, zt);
      const auto fprj = fprj_interp(lt, zt);

      const auto exptau =
        std::exp(tau * (2.0 * mu + tau * sigma * sigma - 2.0 * lc) / 2.0);
      const auto invsqrt_sigma = std::sqrt(2.0) * sigma;
      return (1.0 - fmsk) * (1.0 - fprj) * y3_cluster::invsqrt2pi() / sigma *
               std::exp(-(lc - mu) * (lc - mu) / (2.0 * sigma * sigma)) +
             0.5 * ((1.0 - fmsk) * fprj * tau + fmsk * fprj / lt) * exptau *
               std::erfc((mu + tau * sigma * sigma - lc) / invsqrt_sigma) +
             0.5 * fmsk / lt *
               (std::erfc((lc - mu) / invsqrt_sigma) -
                std::erfc((lc + lt - mu) / invsqrt_sigma)) -
             0.5 * fmsk / lt * fprj *
               (std::exp(-tau * lt) * exptau *
                std::erfc((mu + tau * sigma * sigma - lc - lt) /
                          invsqrt_sigma));
    }
  };
}

#endif
