#ifndef Y3_CLUSTER_OMEGA_Z_SDSS_HH
#define Y3_CLUSTER_OMEGA_Z_T_HH

#include <array>
#include <cmath>

namespace y3_cluster {
  struct  OMEGA_Z_SDSS {

    double
    operator()(double zt) const
    {
      std::array<double, 12> constexpr poly_coeff_vol = {-1.14293122E05,
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
      double omega_z = 0.0;
      double constexpr zpivot = 0.2;

      for (std::size_t i = 0; i < poly_coeff_vol.size(); ++i) {
        omega_z = omega_z +
                  poly_coeff_vol[i] * std::pow(zt - zpivot,  poly_coeff_vol.size() - i - 1.);
      }
      return omega_z;
    }
  };
}

#endif