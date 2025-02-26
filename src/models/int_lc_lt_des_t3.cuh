#ifndef Y3_CLUSTER_INT_LC_LT_DES_T3_CUH
#define Y3_CLUSTER_INT_LC_LT_DES_T3_CUH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "common/cuda/Interp2D.cuh"
#include "utils/make_interp_2d.cuh"
#include "utils/primitives.hh"

#include <algorithm>

namespace y3_cuda {
  class INT_LC_LT_DES_T3 {
  private:
    quad::Interp2D bin1;
    quad::Interp2D bin2;
    quad::Interp2D bin3;
    quad::Interp2D bin4;

  public:
    using doubles = std::vector<double>;

    explicit INT_LC_LT_DES_T3(cosmosis::DataBlock& sample)
      : bin1(make_Interp2D(sample,
                           "grabLoLt",
                           "lt",
                           "zt",
                           "bin1")),
        bin2(make_Interp2D(sample,
                           "grabLoLt",
                           "lt",
                           "zt",
                           "bin2")),
        bin3(make_Interp2D(sample,
                           "grabLoLt",
                           "lt",
                           "zt",
                           "bin3")),                           
        bin4(make_Interp2D(sample,
                           "grabLoLt",
                           "lt",
                           "zt",
                           "bin4"))
    {}

    __device__ __host__ double
    operator()(double lc, double lt, double zt) const
    {
      double val = 0;
      if ((lc >= 20) & (lc < 30)) {
        val = bin1(lt, zt);
      } else if ((lc >= 30) & (lc < 45)) {
        val = bin2(lt, zt);
      } else if ((lc >= 45) & (lc < 60)) {
        val = bin3(lt, zt);
      } else {
        val = bin4.clamp(lt, zt);
      }
      return val;
    }
  };
}

#endif
