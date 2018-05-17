#ifndef Y3_CLUSTER_T_MIS_T_HH
#define Y3_CLUSTER_T_MIS_T_HH

namespace y3_cluster {
  struct T_MIS_t {
    double
    operator()(double, double, double) const
    {
      return 1.0;
    }
  };
}

#endif