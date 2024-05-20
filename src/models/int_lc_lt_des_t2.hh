#ifndef Y3_CLUSTER_INT_LC_LT_DES_T2_HH
#define Y3_CLUSTER_INT_LC_LT_DES_T2_HH

#include "cosmosis/datablock/datablock.hh"
#include "utils/read_vector.hh"
#include "fmt/core.h"

using y3_cluster::Interp2D;

namespace y3_cluster {
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
    
    double
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
  private:
    // CUDA/C++ does not support static data members of classes or structs
    Interp2D const lambda0_interp;
    Interp2D const lambda1_interp;
    Interp2D const lambda2_interp;
    Interp2D const lambda3_interp;
  };
}
#endif
