/* Here I use Omega_m=0.3, Omega_l=0.7 , in class EZ */
/* We would need redshift of the cluster as an input?! */
#include <cmath>

class nfw_xi {
public:
  explicit nfw_xi(double c, double om) : _c(c), _om(om) {}

  explicit nfw_xi(cosmosis::DataBlock& sample)
    : _c(sample.view<double>("cluster_mass_profile", "concentration"))
    , _om(sample.view<double>("cosmological_parameters", "omega_m"))
  {}

  double
  operator()(double r, double M) const /* M in h^-1 M_solar, M_{200} here */
  {
    double const delta_c =
      200. * _c * _c * _c / 3. / (std::log(1. + _c) - _c / (1. + _c));

    double const rho_crit = 2.77536627E11 * _om;
    // this needs to be rho_m rather than rho_crit which is rho_crit*om

    double const r_200 = std::cbrt(3. * M / (800. * M_PI * rho_crit));
    double const r_s = r_200 / _c;
    double const r_ratio = r / r_s;
    return delta_c / (r_ratio * (1 + r_ratio) * (1 + r_ratio)) - 1.0;
  }

private:
  double _c;
  double _om;
};
