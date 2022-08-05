#ifndef Y3_CLUSTER_INT_LC_LT_DES_T_CUH
#define Y3_CLUSTER_INT_LC_LT_DES_T_CUH

#include "cosmosis/datablock/datablock.hh"
#include "cuda/pagani/quad/GPUquad/Interp2D.cuh"
//#include "models/interpolation_tables.cuh"

namespace y3_cuda {
  struct INT_LC_LT_DES_t {

     static quad::Interp2D const lambda0_interp;
     static quad::Interp2D const lambda1_interp;
     static quad::Interp2D const lambda2_interp;
     static quad::Interp2D const lambda3_interp;
    
    explicit INT_LC_LT_DES_t(const cosmosis::DataBlock&) {}
    INT_LC_LT_DES_t()/* :lambda0_interp(lt_bins, zt_bins, lambda0_arr), lambda1_interp(lt_bins, zt_bins, lambda1_arr), lambda2_interp(lt_bins, zt_bins, lambda2_arr), lambda3_interp(lt_bins, zt_bins, lambda3_arr)*/ {}

    __device__ __host__ double
    operator()(double lc, double lt, double zt) const
    {
      double val = 0;
      if ((lc >= 20) & (lc < 30)) {
        val = lambda0_interp(lt, zt);
      } else if ((lc >= 30) & (lc < 45)) {
        val = lambda1_interp(lt, zt);
      } else if ((lc >= 45) & (lc < 60)) {
        val = lambda2_interp(lt, zt);
      } else {
        val = lambda3_interp.clamp(lt, zt);
      }
      return val;
    }
  };
}

#endif
