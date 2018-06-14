#ifndef Y3_CLUSTER_DEL_SIG_MIS_T_HH
#define Y3_CLUSTER_DEL_SIG_MIS_T_HH

namespace y3_cluster {
  struct DEL_SIG_MIS_t {
    double
    operator()(double, double, double) const
    {
      /* eq. (49) */
      // TODO
      return 1.0;
    }
  };
}

#endif
