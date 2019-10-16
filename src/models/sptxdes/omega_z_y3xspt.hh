#ifndef Y3_CLUSTER_OMEGA_Z_Y3XSPT_HH
#define Y3_CLUSTER_OMEGA_Z_Y3XSPT_HH

#include "cosmosis/datablock/datablock.hh"
#include "utils/polynomial.hh"

#include <array>
#include <cmath>

namespace y3_cluster {
  struct OMEGA_Z_Y3XSPT {
  public:
    OMEGA_Z_Y3XSPT() {}
    OMEGA_Z_Y3XSPT(cosmosis::DataBlock&) {}

    double
    operator()(double zt) const
    {
      // These fits will need to be computed by Eli (+Lindsey)
      // For now, use constant area = 2500 deg^2

      // Make a variable for zt for now to avoid unused variable warning
      double const temp_zt = zt;

      // 2500 deg^2 in rad^2
      return 0.7615435494667714;

      // ====================================================================
      // TBD: I need these fits to be calculated for Y3 x SPT-SZ
      // TBD: Also for Y3 x SPT-SZ+SPTpol-ECS
      // TBD: A+ naming SDSS_fit for DES =P
      // // Example: y3_cluster::polynomial<3> f{3,-1,2.5}
      // // yields the function f(x) = 3x^2 - x + 2.5
      // static const y3_cluster::polynomial<6>
      // SDSS_fit{{0.0,0.0,0.0,-0.00262353,
      //                                                  0.01940118,0.45133063}};
      // static const y3_cluster::polynomial<6>
      // SDSS_fit2{{1.33647377e+4,1.35291046e+3,
      //                                                   -1.26204891e+2,-2.83454918e+1,
      //                                                   -2.26465905,3.84958753e-1}};
      // static const y3_cluster::polynomial<6>
      // SDSS_fit3{{0,0,-1.88101967,4.8071839,
      //                                                   -4.11424324,1.18196785}};
      //
      // // Returns effective survey area in rad^2
      // if (zt < 0.504) {return SDSS_fit(zt);}
      // else if (zt < 0.7) {return SDSS_fit2(zt-0.6);}
      // else {return SDSS_fit3(zt);}
    }
  };
}

#endif
