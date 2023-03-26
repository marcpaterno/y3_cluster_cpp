#ifndef Y3_CLUSTER_CPP_MOR_T_SZ_CUH
#define Y3_CLUSTER_CPP_MOR_T_SZ_CUH

#include "cosmosis/datablock/datablock.hh"
#include "utils/mz_sz_power_law.cuh"
#include "utils/primitives.cuh"
#include <cmath>

namespace y3_cuda {

  class MOR_t_sz {
  public:
    MOR_t_sz(mz_sz_power_law lambda, mz_sz_power_law zeta, double sigma1, double sigma2, double corr12)
      : _lambda(lambda), _zeta(zeta), _sigma_intr1(sigma1), _sigma_intr2(sigma2), _corr12(corr12)
    {}
    explicit MOR_t_sz(cosmosis::DataBlock& sample)
      : _lambda(sample.view<double>("scalingRelation", "A1"),
                sample.view<double>("scalingRelation", "B1"),
                sample.view<double>("scalingRelation", "C1"),
                sample.view<double>("cosmological_parameters", "omega_M"))

      , _zeta(sample.view<double>("scalingRelation", "A2"),
              sample.view<double>("scalingRelation", "B2"),
              sample.view<double>("scalingRelation", "C2"),
              sample.view<double>("cosmological_parameters", "omega_M"))

      , _sigma_intr1(sample.view<double>("scalingRelation", "sigma1"))
      , _sigma_intr2(sample.view<double>("scalingRelation", "sigma2"))
      , _corr12(sample.view<double>("scalingRelation", "corr12"))
    {}

    __device__ __host__ double
    operator()(double lt, double szt, double lnM, double zt) const
    {
      /* eq. 2-3 de Grandis et al. 2021 */
      double const ltm = _lambda(lnM, zt);
      double const sztm = _zeta(lnM, zt);

      double const _sigma1 =
        max(sqrt(_sigma_intr1 * _sigma_intr1  + 1/ltm), 4.0);
      double const _sigma2 = max(_sigma_intr2, 0.5);
      double const _corr = min(max(_corr12, 1.),-1.);

      double const x1 = log(lt) - log(ltm);
      double const x2 = log(szt) - log(sztm);
      double const res = y3_cuda::gaussian2d(x1, x2, 0.0, 0.0, _sigma1, _sigma2, _corr);
      return res;
    }

  private:
    mz_sz_power_law _lambda;
    mz_sz_power_law _zeta;
    double _sigma_intr1;
    double _sigma_intr2;
    double _corr12;
  };
}

#endif
