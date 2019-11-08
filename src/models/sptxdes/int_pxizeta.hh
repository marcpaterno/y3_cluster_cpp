#ifndef Y3_CLUSTER_INT_PXIZETA_HH
#define Y3_CLUSTER_INT_PXIZETA_HH

#include "utils/datablock_reader.hh"
#include <cmath>

namespace y3_cluster {

  class INT_PXIZETA {
  public:
    explicit INT_PXIZETA() {}

    double
    operator()(double zeta) const
    {
      // The integral over xi in exp(-(xi-meanzeta)**2/2)/sqrt(2pi)
      // meanzeta = sqrt(zeta^2 + 3)
      // I(zeta) = 1/2 erfc(-sqrt(zeta^2+3)/sqrt(2))
      return std::erfc(-std::sqrt((zeta * zeta + 3) / 2.0)) / 2.0;
    }
  };
}

#endif
