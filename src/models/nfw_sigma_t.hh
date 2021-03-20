/* Here I use Omega_m=0.3, Omega_l=0.7 , in class EZ */
/* We would need redshift of the cluster as an input?! */
#include <iostream>
#include <math.h>
namespace y3_cluster {
  class nfw_sigma {
  public:
    explicit nfw_sigma(double c, double om) : _c(c), _om(om) {}
    using doubles = std::vector<double>;
    explicit nfw_sigma(cosmosis::DataBlock& sample)
      : _c(sample.view<double>("cluster_mass_profile", "concentration"))
      , _om(sample.view<double>("cosmological_parameters", "omega_m"))
    {}

    double
    operator()(double r, double M) const /* M in h^-1 M_solar, M_{200} here */
    {

      // Wright and Brainerd 1999
      double const delta_c =
        200. * _c * _c * _c / 3. * (std::log(1. + _c) - _c / (1. + _c));

      double const rho_crit = 2.77536627E11 * _om;
      /* EZ*EZ would be a little bit slower than direct definition */

      double const r_200 = std::cbrt(3. * M / (800. * M_PI * rho_crit));

      double const r_s = r_200 / _c;

      double const r_ratio = r / r_s;

      double const coeff = r_s * delta_c * rho_crit / 1E12;

      if (r_ratio < 1.0) {
        return  2.0 * coeff / (r_ratio * r_ratio - 1) *
              (1 - 2.0 / std::sqrt(1 - r_ratio * r_ratio) *
                     std::atanh(std::sqrt((1. - r_ratio) / (1. + r_ratio))));
      }
      if (r_ratio == 1.0) {
        return  2.0 * coeff / 3.0;
      }

      return 2.0 * coeff / (r_ratio * r_ratio - 1) *
              (1 - 2.0 / std::sqrt(r_ratio * r_ratio - 1) *
                     std::atan(std::sqrt((r_ratio - 1) / (1. + r_ratio))));
    }

  private:
    double _c;
    double _om;
  };
}
