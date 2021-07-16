#ifndef EZ_CUH
#define EZ_CUH

#include "models/ez_sq.cuh"

class ez {
public:
  ez() = default;

  ez(double omega_m, double omega_l, double omega_k)
    : _ezsq(omega_m, omega_l, omega_k)
  {}

  __device__ double
  operator()(double z) const
  {
    auto const sqr = _ezsq(z);
    return sqrt(sqr);
  }

private:
  ez_sq _ezsq;
};

#endif
