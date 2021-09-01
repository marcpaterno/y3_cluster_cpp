#ifndef Y3_CLUSTER_CPP_MOR_T_CUH
#define Y3_CLUSTER_CPP_MOR_T_CUH

#include "cosmosis/datablock/datablock.hh"
#include "utils/mz_power_law.cuh"
#include "utils/primitives.cuh"
#include <cmath>

class MOR_t {
public:
  MOR_t(mz_power_law lambda, double sigma_i, double alpha)
    : _lambda(lambda), _sigma_intr(sigma_i), _alpha(alpha)
  {}

  explicit MOR_t(cosmosis::DataBlock& sample)
    : _lambda(sample.view<double>("cluster_abundance", "mor_A"),
              sample.view<double>("cluster_abundance", "mor_B"),
              sample.view<double>("cluster_abundance", "mor_epsilon"))
    , _sigma_intr(sample.view<double>("cluster_abundance", "mor_sigma"))
    , _alpha(sample.view<double>("cluster_abundance", "mor_alpha"))
  {}

  __device__ __host__ double
  operator()(double lt, double lnM, double zt) const
  {
    /* eq. (34) */
    double const ltm = _lambda(lnM, zt);
    double const _sigma =
      max(sqrt(_sigma_intr * _sigma_intr * ltm * ltm + ltm), 4.0);
    double const x = lt - ltm;
    double const erfarg = -1.0 * _alpha * (x) / (sqrt(2.) * _sigma);
    double const erfterm = erfc(erfarg);
    return y3_cuda::gaussian(x, 0.0, _sigma) * erfterm;
  }

private:
  mz_power_law _lambda;
  double _sigma_intr;
  double _alpha;
};

#endif
