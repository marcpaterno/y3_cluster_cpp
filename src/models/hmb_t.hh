#ifndef Y3_CLUSTER_CPP_BMZ_T_HH
#define Y3_CLUSTER_CPP_BMZ_T_HH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "cosmosis/datablock/datablock.hh"
#include "utils/interp_2d.hh"
#include "utils/make_interp_2d.hh"
#include "utils/read_vector.hh"

#include <memory>

namespace y3_cluster {
  // Bias for the Halo Mass Function (Halo Mass Bias)
  class HMB_t {
    Interp2D _tinker_bias;

  public:
    HMB_t()
      : _tinker_bias(read_vector("tinker_bias_function/ln_mass_h.txt"),
                     read_vector("tinker_bias_function/z.txt"),
                     read_vector("tinker_bias_function/bias.txt"))
    {}

    explicit HMB_t(cosmosis::DataBlock& sample)
      : _tinker_bias(make_Interp2D(sample, "tinker_bias_function", "ln_mass_h","z", "bias"))
    {}

    double
    operator()(double lnM, double zt) const
    {
      return _tinker_bias(lnM, zt);
    }
  };
}

#endif
