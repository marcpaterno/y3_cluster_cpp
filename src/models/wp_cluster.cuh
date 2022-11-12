#ifndef Y3_CLUSTER_XI_SUM_CUH
#define Y3_CLUSTER_XI_SUM_CUH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "cuda/pagani/quad/GPUquad/Interp2D.cuh"
#include "utils/make_interp_2d.cuh"
#include "utils/primitives.cuh"

namespace y3_cuda {
  class WP_CLUSTER {
  private:
    quad::Interp2D _sigma2;
    quad::Interp2D _bias;

  public:
    WP_CLUSTER(quad::Interp2D const& sigma2,
            quad::Interp2D const& bias)
      :  _sigma2(sigma2), _bias(bias)
    {}

    using doubles = std::vector<double>;

    explicit WP_CLUSTER(cosmosis::DataBlock& sample)
      , _sigma2(make_Interp2D(sample,
                              "deltasigma",
                              "r_sigma_deltasigma",
                              "matter_power_lin",
                              "z",
                              "deltasigma",
                              "sigma_2"))
      , _bias(make_Interp2D(sample,
                            "matter_power_lin",
                            "z",
                            "deltasigma",
                            "lnM",
                            "deltasigma",
                            "bias"))
    {}

    __device__ __host__ double
    operator()(double r, double lnM, double zt) const
    // The cluster cluster correlation function is the FT(PS)*b^2
    // so is the same as the 2halo term in sigma profile with another bias term
    /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
    {
      double const sig_2 = _bias.clamp(zt, lnM) * _bias.clamp(zt, lnM) * _sigma2.clamp(r, zt);
      return sig_2;
    }
  };
}

#endif
