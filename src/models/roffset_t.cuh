#ifndef ROFFSET_T_CUH
#define ROFFSET_T_CUH

#include "cosmosis/datablock/datablock.hh"

class ROFFSET_t {
public:
  ROFFSET_t() = default;

  ROFFSET_t(cosmosis::DataBlock& sample)
    : _tau(sample.view<double>("cluster_abundance", "roffset_tau"))
  {}

  explicit ROFFSET_t(double tau) : _tau(tau) {}

  __device__ double
  operator()(double x) const
  {
    // eq. 36
    return x / _tau / _tau * exp(-x / _tau);
  }

private:
  double _tau = 0.0;
};

#endif
