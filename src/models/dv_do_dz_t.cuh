#ifndef Y3_CLUSTER_DV_DO_DZ_T_CUH
#define Y3_CLUSTER_DV_DO_DZ_T_CUH

#include "cosmosis/datablock/datablock.hh"
#include "cudaPagani/quad/GPUquad/Interp1D.cuh"
#include "ez.cuh"
#include "utils/datablock_reader.hh"
#include "utils/make_interp_1d.cuh"

namespace y3_cuda {
  class DV_DO_DZ_t {
  public:
    explicit DV_DO_DZ_t(cosmosis::DataBlock& sample)
      : _da(make_Interp1D(sample, "distances", "z", "d_a"))
      , _ezt(EZ(sample))
      , _h(get_datablock<double>(sample, "cosmological_parameters", "h0"))
    {}

    __device__ __host__ double
    operator()(double zt) const
    {
      double const da_z = _da(zt); // da_z needs to be in Mpc
      // Units: (Mpc/h)^3
      // 2997.92 is Hubble distance, c/H_0
      return 2997.92 * (1.0 + zt) * (1.0 + zt) * da_z * _h * da_z * _h /
             _ezt(zt);
    }

  private:
    quad::Interp1D _da;
    EZ _ezt;
    double _h;
  };
}

#endif
