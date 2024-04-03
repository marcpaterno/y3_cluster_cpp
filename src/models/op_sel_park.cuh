#ifndef Y3_CLUSTER_OP_SEL_PARK_CUH
#define Y3_CLUSTER_OP_SEL_PARK_CUH

#include "cosmosis/datablock/datablock.hh"
#include "common/cuda/Interp2D.cuh"
#include <vector>

namespace y3_cuda {
  // Function to calculate the Park et al. (2023) pi function 
  class OP_SEL_PARK {
  public:
    // Initialize directly with scalar values for PI0, R0, and c
    OP_SEL_PARK(double PI0_scalar, double R0_scalar, double c_scalar)
      : _PI0_val(PI0_scalar), _R0_val(R0_scalar), _c(c_scalar),
        _PI0(PI0_scalar, 0.0, 0.0), // Direct initialization with A, B, and C
        _R0(R0_scalar, 0.0, 0.0)    // Assuming B=0.0 and C=0.0 as placeholders
    {}

    // Initialize using cosmosis::DataBlock for data-driven computation
    explicit OP_SEL_PARK(cosmosis::DataBlock& sample)
      : _PI0_val(sample.view<double>("park_model", "PI0_a")),
        _R0_val(sample.view<double>("park_model", "R0_a")),
        _c(sample.view<double>("park_model", "c")),
        _PI0(sample.view<double>("park_model", "PI0_a"), 
             sample.view<double>("park_model", "PI0_b"), 
             sample.view<double>("park_model", "PI0_c")),
        _R0(sample.view<double>("park_model", "R0_a"), 
            sample.view<double>("park_model", "R0_b"), 
            sample.view<double>("park_model", "R0_c"))
    {}

    __device__ __host__ double operator()(double r, double x, double z) const {
      double const pi0 = _PI0(x,z);
      double const r0 = _R0(x,z);
      return park_pi_func(r, r0, pi0, _c);
    }

    __device__ __host__ double operator()(double r) const {
      return park_pi_func(r, _R0_val, _PI0_val, _c);
    }

  private:
    double const _PI0_val;
    double const _R0_val;
    y3_cuda::par_power_law _PI0;
    y3_cuda::par_power_law _R0;
    double const _c;
  };
}

#endif
