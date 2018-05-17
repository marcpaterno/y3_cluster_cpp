#ifndef Y3_CLUSTER_T_CEN_T_HH
#define Y3_CLUSTER_T_CEN_T_HH

namespace y3_cluster {

  struct T_CEN_t {
    double
    operator()(double, double) const
    {
      return 1.0;
    }
  };
}
#endif