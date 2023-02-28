#ifndef Y3_CLUSTER_INT_LC_LT_DES_T_CUH
#define Y3_CLUSTER_INT_LC_LT_DES_T_CUH

#include "cosmosis/datablock/datablock.hh"
#include "common/cuda/Interp2D.cuh"

namespace y3_cuda {
  struct INT_LC_LT_DES_t {

    // CUDA/C++ does not support static data members of classes or structs
    quad::Interp2D lambda0_interp;
    quad::Interp2D lambda1_interp;
    quad::Interp2D lambda2_interp;
    quad::Interp2D lambda3_interp;

    // The default constructor is implemented in the .cu file.
    INT_LC_LT_DES_t();

    // The constructor from a datablock is required, but does the same thing
    // as the default constructor. So we delegate the work to the default
    // constructor.
    explicit INT_LC_LT_DES_t(cosmosis::DataBlock const&) : INT_LC_LT_DES_t() { }

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
