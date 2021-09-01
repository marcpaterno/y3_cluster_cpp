#ifndef EZ_CUH
#define EZ_CUH

#include "models/ez_sq.cuh"

class EZ {
public:

  EZ(double omega_m, double omega_l, double omega_k)
    : _ezsq(omega_m, omega_l, omega_k)
  {}

  explicit EZ(cosmosis::DataBlock& sample)
    : EZ(sample.view<double>("cosmological_parameters", "omega_m"),
        sample.view<double>("cosmological_parameters", "omega_lambda"),
        sample.view<double>("cosmological_parameters", "omega_k"))
  {}

  __device__ __host__ double
  operator()(double z) const
  {
    auto const sqr = _ezsq(z);
    return sqrt(sqr);
  }

private:
  EZ_sq _ezsq;
};

#endif
