#ifndef Y3_CLUSTER_CPP_PXIZETA_T_HH
#define Y3_CLUSTER_CPP_PXIZETA_T_HH

#include "utils/primitives.hh"

#include <cmath>

namespace y3_cluster {

  class PXiZeta_t {
  public:
    PXiZeta_t() {}

    double
    operator()(double xi, double zeta) const
    {
      // Observational uncertainty in SZ significance P(\xi|\zeta)
      // Eq. 6 of Bleem et al. 2019 (IR version)
      double const mean = std::sqrt(zeta * zeta + 3.0);
      return y3_cluster::gaussian(xi, mean, 1.0);
    }
  };

}

#endif
