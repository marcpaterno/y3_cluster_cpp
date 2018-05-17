#ifndef Y3_CLUSTER_A_MIS_T_HH
#define Y3_CLUSTER_A_MIS_T_HH

namespace y3_cluster {
  struct A_MIS_t {
    double
    operator()(double, double, double, double, double) const
    {
      return 1.0;
    }
  };
}

#endif