#ifndef Y3_CLUSTER_INT_LC_LT_DES_T2_CUH
#define Y3_CLUSTER_INT_LC_LT_DES_T2_CUH

#include "cosmosis/datablock/datablock.hh"
// #include "common/cuda/Interp2D.cuh"
#include "utils/cuda_interp_2d.cuh"
#include "utils/read_vector.hh"
#include "fmt/core.h"

namespace y3_cuda {
  // Define the interpolator values for the LC_LT_DES_t model
  // zt_bins, lt_bins, lambda{bin}_arr
    static inline std::string
    lambda_file(std::string const& BIN_ID)
    {
      return fmt::format("int_lc_lt_des_t_lbd_{}.txt",
                         BIN_ID);
    }

  class INT_LC_LT_DES_t2 {
  public:
    // CUDA/C++ does not support static data members of classes or structs
    gpu_support::Interp2D lambda0_interp;
    gpu_support::Interp2D lambda1_interp;
    gpu_support::Interp2D lambda2_interp;
    gpu_support::Interp2D lambda3_interp;

    INT_LC_LT_DES_t2()
      : lambda0_interp(read_vector("int_lc_lt_des_t_lt_bins.txt"),
                       read_vector("int_lc_lt_des_t_zt_bins.txt"),
                       read_vector(lambda_file("bin0"))),
        lambda1_interp(read_vector("int_lc_lt_des_t_lt_bins.txt"),
                       read_vector("int_lc_lt_des_t_zt_bins.txt"),
                       read_vector(lambda_file("bin1"))),
        lambda2_interp(read_vector("int_lc_lt_des_t_lt_bins.txt"),
                       read_vector("int_lc_lt_des_t_zt_bins.txt"),
                       read_vector(lambda_file("bin2"))),
        lambda3_interp(read_vector("int_lc_lt_des_t_lt_bins.txt"),
                       read_vector("int_lc_lt_des_t_zt_bins.txt"),
                       read_vector(lambda_file("bin3")))
                       { }
    
    explicit INT_LC_LT_DES_t2(cosmosis::DataBlock const&) : INT_LC_LT_DES_t2() { }

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
  // private:
  //   // CUDA/C++ does not support static data members of classes or structs
  //   gpu_support::Interp2D lambda0_interp;
  //   gpu_support::Interp2D lambda1_interp;
  //   gpu_support::Interp2D lambda2_interp;
  //   gpu_support::Interp2D lambda3_interp;
  };
}

#endif
