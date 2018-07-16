#ifndef Y3_CLUSTER_CPP_EZ_HH
#define Y3_CLUSTER_CPP_EZ_HH

#include "ez_sq.hh"
#include <cmath>

namespace y3_cluster {
  class EZ {
  public:
    EZ(double omega_m, double omega_l, double omega_k)
      : _ezsq(omega_m, omega_l, omega_k)
    {}

    double
    operator()(double z) const
    {
      auto const sqr = _ezsq(z);
      return std::sqrt(sqr);
    }

  private:
    EZ_sq _ezsq;
  };
}

#endif
