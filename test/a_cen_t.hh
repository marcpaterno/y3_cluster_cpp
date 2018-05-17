#ifndef Y3_CLUSTER_A_CEN_T_HH
#define Y3_CLUSTER_A_CEN_T_HH

namespace y3_cluster {
  struct A_CEN_t {
    double
    operator()(double, double, double, double) const
    {
      return 1.0;
    }
  };
}

#endif