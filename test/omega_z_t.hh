#ifndef Y3_CLUSTER_OMEGA_Z_T_HH
#define Y3_CLUSTER_OMEGA_Z_T_HH

#include <cmath>

namespace y3_cluster {
  class OMEGA_Z_t {
  public:
    double
    operator()(double zt) const
    {
      double poly_coeff_vol[12] = {-1.14293122E05,
                                   5.96846869E04,
                                   9.24239180E03,
                                   -2.23118813E03,
                                   -4.52580713E03,
                                   1.18404878E03,
                                   1.27951911E02,
                                   -5.05716847E01,
                                   1.01744577E00,
                                   -3.11253383E-01,
                                   5.48481084E-03,
                                   3.12629987E00};
      int poly_deg = 12;
      double omega_z = 0.0;
      double zpivot = 0.2;

      for (int i = 0; i < 12; i++) {
        omega_z = omega_z +
                  poly_coeff_vol[i] * std::pow(zt - zpivot, poly_deg - i - 1.);
      }
      return omega_z;
    }
  };
}

#endif