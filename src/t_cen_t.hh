#ifndef Y3_CLUSTER_T_CEN_T_HH
#define Y3_CLUSTER_T_CEN_T_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"

namespace y3_cluster {

  class T_CEN_t {
    double const ln_r0;
    double const gamma;
  public:
    // hard-coded values are best-fits from triaxiality paper
    T_CEN_t()
      : ln_r0(0.632)
      , gamma(1.634)
    {}
    explicit T_CEN_t(cosmosis::DataBlock&)
      : T_CEN_t()
    {}

    double
    operator()(double r, double lnM [[maybe_unused]]) const
    {
      // FIXME: What should units of r be before log? Mpc, Mpc/h?
      // Consult triaxiality paper
      double const ln_r = std::log(r),
                   ln_r_diff = ln_r0 - ln_r;
      return 1 - 1.0 / ((ln_r_diff * ln_r_diff) + gamma);
    }
  };
}

#endif
