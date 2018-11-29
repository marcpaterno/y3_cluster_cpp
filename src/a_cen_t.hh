#ifndef Y3_CLUSTER_A_CEN_T_HH
#define Y3_CLUSTER_A_CEN_T_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "polynomial.hh"

namespace y3_cluster {

  class A_CEN_t {
  public:
      A_CEN_t() {}
      explicit A_CEN_t(cosmosis::DataBlock&) {}

      double
      operator()(double mu) const
      {
        // Coefficients are best-fit vals from
        // "Effects of Cluster Triaxiality on DES Stacked Weak Lensing Analysis"
        // Zhang, Wu, and Zhang, (in preparation)
        static constexpr polynomial<4> acen{{{0.637, -0.191, 0.182, -0.191}}};
        return acen(mu);
      }
  };

}

#endif
