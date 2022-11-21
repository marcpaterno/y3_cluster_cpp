#ifndef Y3_CLUSTER_PK_CLUSTER_DAMP_CUH
#define Y3_CLUSTER_PK_CLUSTER_DAMP_CUH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "cuda/pagani/quad/GPUquad/Interp2D.cuh"
#include "utils/make_interp_2d.cuh"
#include "utils/primitives.cuh"

namespace y3_cuda {
  class PK_CLUSTER_DAMP {
  private:
    quad::Interp2D _pk_damp;
    quad::Interp2D _bias;

  public:
    size_t
    get_device_mem_footprint()
    {
      size_t size = 0;
      size += _pk_damp.get_device_mem_footprint();
      size += _bias.get_device_mem_footprint();
      return size;
    }

    PK_CLUSTER_DAMP(quad::Interp2D const& pk_damp,
            quad::Interp2D const& bias)
      :  _pk_damp(pk_damp), _bias(bias)
    {}

    using doubles = std::vector<double>;

    explicit PK_CLUSTER_DAMP(cosmosis::DataBlock& sample)
      : _pk_damp(make_Interp2D(sample,
                              "correlationFunction",
                              "k",
                              "z",
                              "Damped_Pk_hh"))
      , _bias(make_Interp2D(sample,
                              "correlationFunction",
                              "z",
                              "lnM",
                              "bias"))
    {}

    __device__ __host__ double
    operator()(double k, double lnM, double zt) const
    // The 3d cluster cluster power spectrum Pk_hh*b^2
    // With a damping term due to the cluster photo-z redshift error
    // so is the matter-matter power spectrum ?plus the NFW power spectrum
    {
      double const pk_damp = _bias.clamp(zt, lnM) * _bias.clamp(zt, lnM) * _pk_damp.clamp(k, zt);
      return pk_damp;
    }
  };
}

#endif
