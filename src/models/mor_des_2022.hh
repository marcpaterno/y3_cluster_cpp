#ifndef Y3_CLUSTER_CPP_MOR_DES_2022_HH
#define Y3_CLUSTER_CPP_MOR_DES_2022_HH

#include "cosmosis/datablock/datablock.hh"
#include "utils/interp_2d.hh"
#include "utils/mz_power_law.hh"
#include "utils/primitives.hh"

#include <cmath>
#include <string>

namespace y3_cluster {

  class MOR_DES_2022 {

    double _A;
    double _B;
    double _C;
    double _sigma_intr;
    double _epsilon;
    double _z_pivot;

  public:
    MOR_DES_2022(double A,
              double B,
              double C,
              double sigma_i,
              double epsilon,
              double z_pivot)
      : _A(A)
      , _B(B)
      , _C(C)
      , _sigma_intr(sigma_i)
      , _epsilon(epsilon)
      , _z_pivot(z_pivot)
    {}

    explicit MOR_DES_2022(cosmosis::DataBlock& sample)
      : _A(sample.view<double>("cluster_abundance", "mor_A"))
      , _B(sample.view<double>("cluster_abundance", "mor_B"))
      , _C(sample.view<double>("cluster_abundance", "mor_alpha"))
      , _sigma_intr(sample.view<double>("cluster_abundance", "mor_sigma"))
      , _epsilon(sample.view<double>("cluster_abundance", "mor_epsilon"))
      , _z_pivot(sample.view<double>("cluster_abundance", "z_mor_pivot"))
    {}

    double
    operator()(double lt, double lnM, double zt) const
    {

      double const l_sat_M = pow((std::exp(lnM) - _A) / (_B - _A), _C) *
                         pow((1.0 + zt) / (1.0 + _z_pivot), _epsilon);
      //l_sat_M =np.exp(LnLambda_M(mass,z_reds))-1.0 // Expected number of satellite galaxies
      double const shift_x = (l_sat_M*_sigma_intr)*(l_sat_M*_sigma_intr); // Gaussian error
      double const std = std::sqrt(l_sat_M+shift_x);  // sum in quadrature poisson and gauss error
      double const ln_gamma_fun = std::lgamma(lt+shift_x);    // ln Gamma function

      return std::exp(-std*std+(lt + shift_x-1.) * std::log(std*std)-ln_gamma_fun);

    }
  };
}

#endif
