#ifndef Y3_CLUSTER_CPP_ROFFSET_T_HH
#define Y3_CLUSTER_CPP_ROFFSET_T_HH

#include <cmath>

namespace y3_cluster {

  class ROFFSET_t {
  public:
    explicit ROFFSET_t(double tau) : _tau(tau) {}

    double
    operator()(double x) const
    {
      // eq. 36
      return x / _tau / _tau * std::exp(-x / _tau);
    }

  private:
    double _tau;
  };
}

#endif
