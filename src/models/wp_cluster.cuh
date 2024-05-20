#ifndef Y3_CLUSTER_WP_SUM_CUH
#define Y3_CLUSTER_WP_SUM_CUH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "common/cuda/Interp2D.cuh"
#include "utils/make_interp_2d.cuh"
#include "utils/primitives.cuh"

// Computes the Average 3D 2PCF: \xi(r)
// Wp is obtained after by projecting \xi(r) onto the plane of the sky
// TODO 
// CHECK THE MODEL DEFINITION

namespace y3_cuda {
  class WP_CLUSTER {
  private:
    quad::Interp2D _wp1;
    quad::Interp2D _wp2;
    quad::Interp2D _bias;

  public:
    size_t
    get_device_mem_footprint()
    {
      size_t size = 0;
      size += _wp1.get_device_mem_footprint();
      size += _wp2.get_device_mem_footprint();
      size += _bias.get_device_mem_footprint();
      return size;
    }

    WP_CLUSTER(quad::Interp2D const& sigma1,
              quad::Interp2D const& sigma2,
              quad::Interp2D const& bias)
      : _wp1(sigma1), _wp2(sigma2), _bias(bias)
    {}

    using doubles = std::vector<double>;

    explicit WP_CLUSTER(cosmosis::DataBlock& sample)
      : _wp1(make_Interp2D(sample,
                              "haloModel",
                              "Rp",
                              "lnM",
                              "Wp_nfw"))
      , _wp2(make_Interp2D(sample,
                              "haloModel",
                              "Rp",
                              "z",
                              "Wp_hh"))
      , _bias(make_Interp2D(sample,
                            "haloModel",
                            "z",
                            "lnM",
                            "bias"))
    {}

    __device__ __host__ double
    operator()(double r, double lnM, double zt) const
    /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
    // TODOs: Check if it's a double integral on Mass for the bias
    // TODOs: \int \int dM1 dM2 bias(M1) bias(M2) ...
    {
      double const wp_1 = _wp1.clamp(r, lnM);
      double const wp_2 = _bias.clamp(zt, lnM) * _bias.clamp(zt, lnM)  * _wp2.clamp(r, zt);
      double const res = wp_1+wp_2;
      return res;
    }
  };
}

#endif
// Jim's implementation wo the 1h term
// #ifndef Y3_CLUSTER_WP_SUM_CUH
// #define Y3_CLUSTER_WP_SUM_CUH

// #include "cosmosis/datablock/datablock.hh"
// #include "cosmosis/datablock/ndarray.hh"
// #include "cuda/pagani/quad/GPUquad/Interp2D.cuh"
// #include "utils/make_interp_2d.cuh"
// #include "utils/primitives.cuh"

// namespace y3_cuda {
//   class WP_CLUSTER {
//   private:
//     quad::Interp2D _sigma2;
//     quad::Interp2D _bias;

//   public:
//     size_t
//     get_device_mem_footprint()
//     {
//       size_t size = 0;
//       size += _sigma2.get_device_mem_footprint();
//       size += _bias.get_device_mem_footprint();
//       return size;
//     }

//     WP_CLUSTER(quad::Interp2D const& sigma2,
//             quad::Interp2D const& bias)
//       :  _sigma2(sigma2), _bias(bias)
//     {}

//     using doubles = std::vector<double>;

//     explicit WP_CLUSTER(cosmosis::DataBlock& sample)
//       : _sigma2(make_Interp2D(sample,
//                               "deltasigma",
//                               "r_sigma_deltasigma",
//                               "matter_power_lin",
//                               "z",
//                               "deltasigma",
//                               "sigma_2"))
//       , _bias(make_Interp2D(sample,
//                             "matter_power_lin",
//                             "z",
//                             "deltasigma",
//                             "lnM",
//                             "deltasigma",
//                             "bias"))
//     {}

//     __device__ __host__ double
//     operator()(double r, double lnM, double zt) const
//     // The cluster cluster correlation function is the FT(PS)*b^2
//     // so is the same as the 2halo term in sigma profile with another bias term
//     /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
//     {
//       double const sig_2 = _bias.clamp(zt, lnM) * _bias.clamp(zt, lnM) * _sigma2.clamp(r, zt);
//       return sig_2;
//     }
//   };
// }

// #endif
