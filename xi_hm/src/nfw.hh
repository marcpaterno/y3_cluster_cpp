/* Here I use Omega_m=0.3, Omega_l=0.7 , in class EZ */
/* We would need redshift of the cluster as an input?! */
#include <math.h>

class nfw_dsigma{
public:
  explicit nfw_dsigma(double c)
    : _c(c)
  {}

  explicit nfw_dsigma(cosmosis::DataBlock& sample)
  { // remember to read this in from the values.ini file
    sample.get_val<double>("nfw_dsigma_params", "c", _c);
  }

  double
  operator()(double r, double M) const /* M in h^-1 M_solar, M_{200} here */
  {
    double delta_c = 200.*_c*_c*_c / 3./(std::log(1.+ _c) - _c/(1.+ _c));

    double rho_crit = 2.77526157E11;
    // this needs to be rho_m rather than rho_crit which is rho_crit*om

    double r_200 = std::pow(3.*M/(800.*M_PI*rho_crit) , 1./3.);
    double r_s = r_200 / _c;
    double r_ratio = r / r_s;
    double xi_nfw = delta_c/(r_ratio*(1+r_ratio)*(1+r_ratio)) - 1.0;

    return xi_nfw; 
 }
private:
  double _c;
};
