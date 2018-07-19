#ifndef Y3_CLUSTER_CPP_BMZ_T_HH
#define Y3_CLUSTER_CPP_BMZ_T_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"

namespace y3_cluster {
  // Bias for the Halo Mass Function
  class BMZ_t {
  public:
    BMZ_t() {}
    explicit BMZ_t(cosmosis::DataBlock&) {}

    double
    operator()(double /* lnM */, double /* zt */) const
    {
      return 1.0;
    }
  };
}

#endif
