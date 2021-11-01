#ifndef Y3_CLUSTER_INT_ZO_ZT_CUH
#define Y3_CLUSTER_INT_ZO_ZT_CUH

#include "models/sigma_photoz_des.cuh"

namespace y3_cuda {
  class INT_ZO_ZT_DES_t {
  public:
    __device__ __host__ double
    operator()(double zomin, double zomax, double zt) const
    {
      double _sigma = _sigma_photoz_des(zt);
      double base = sqrt(2) * _sigma;
      return (erf((zomax - zt) / base) - erf((zomin - zt) / base)) / 2.0;
    }

  private:
    SIGMA_PHOTOZ_DES_t _sigma_photoz_des;
  };
}

#endif
