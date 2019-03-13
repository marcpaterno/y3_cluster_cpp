/* Here I use Omega_m=0.3, Omega_l=0.7 , in class EZ */
/* We would need redshift of the cluster as an input?! */
#include <math.h>

class halo_2{
public:
  explicit  halo_2(double c)
    : _c(c)
  {}

  explicit halo_2(cosmosis::DataBlock& sample)
  {
    sample.get_val<double>("nfw_dsigma_params", "c", _c);
  }

  double
  operator()(double r, double lnM) const /* M in h^-1 M_solar, M_{200} here */
  {
    //EZ ez{0.3 , 0.7, 0.0};
    double ez = 1.0;
    double delta_c = 200.*_c*_c*_c / 3.*(std::log(1.+ _c) - _c/(1.+ _c));

    double rho_crit = 2.77526157E11*ez*ez;
    /* EZ*EZ would be a little bit slower than direct definition */

    double r_200 = std::pow(3.*std::exp(lnM)/(800.*M_PI*rho_crit) , 1./3.);

    double r_s = r_200 / _c;

    double r_ratio = r / r_s;

    double coeff = r_s * delta_c * rho_crit;

    if(r_ratio<1.){
      return coeff * ( 8.* std::atanh(std::sqrt((1.-r_ratio)/(1.+r_ratio)))/
              (r_ratio*r_ratio*std::sqrt(1.-r_ratio*r_ratio))
              +4.*std::log(r_ratio/2.)/(r_ratio*r_ratio)
              -2./(r_ratio*r_ratio-1.)
              +4.*std::atanh(std::sqrt((1.-r_ratio)/(1.+r_ratio)))/
              ((r_ratio*r_ratio-1.)*std::sqrt(1.-r_ratio*r_ratio)) );
    }else if(r_ratio==1.){
      return coeff * (10./3. + 4.*std::log(0.5));
    }else{
      return coeff * ( 8.* std::atan(std::sqrt((r_ratio-1.)/(r_ratio+1.)))/
              (r_ratio*r_ratio*std::sqrt(r_ratio*r_ratio-1.))
              +4.*std::log(r_ratio/2.)/(r_ratio*r_ratio)
              -2./(r_ratio*r_ratio-1.)
              +4.*std::atanh(std::sqrt((r_ratio-1.)/(r_ratio+1.)))/
              (std::pow(r_ratio*r_ratio-1.,1.5)) );
    }
 }
private:
  double _c;
};
