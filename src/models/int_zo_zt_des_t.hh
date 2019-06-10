#ifndef Y3_CLUSTER_INT_ZO_ZT_HH
#define Y3_CLUSTER_INT_ZO_ZT_HH

#include "utils/datablock_reader.hh"
#include <cmath>

namespace y3_cluster {

  class INT_ZO_ZT_DES_t {
  public:
    explicit INT_ZO_ZT_DES_t(){}


    double
    operator()(double zomin, double zomax, double zt) const
    {

      double poly_coeff[]={ -40358.8315,   2798.08304,   9333.80185,  -657.348248,  -840.565610,   46.8506649,   37.8839498,  -0.868811858,  -0.808928182,   0.00890199353,   0.0139811265};
      using std::erf;
      double _sigma=0;
      double z_for_fit=zt;
      if (z_for_fit<0.15) { z_for_fit=0.15; } // no extrapolation outside data range
      if (z_for_fit>0.7)  { z_for_fit=0.7; }// no extrapolation outside data range
	
      z_for_fit = z_for_fit - 0.4;      
      for (int ii=0;ii<10;ii++){      
      _sigma = (poly_coeff[ii]+_sigma)*z_for_fit; 
      }
      _sigma=_sigma+poly_coeff[10];
      double base = std::sqrt(2) * _sigma;
      return (erf((zomax - zt) / base) -
              erf((zomin - zt) / base)) / 2.0;
    }

  private:
  };
}

#endif
