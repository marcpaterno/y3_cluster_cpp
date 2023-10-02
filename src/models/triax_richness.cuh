#ifndef Y3_CLUSTER_TRIAX_RICHNESS_CUH
#define Y3_CLUSTER_TRIAX_RICHNESS_CUH

#include "cosmosis/datablock/datablock.hh"
#include <array>
#include "common/cuda/Interp1D.cuh"

namespace y3_cuda {
  class TRIAX_RICHNESS {
  private:
    quad::Interp1D _lnA;

  public:
    /*The values are taken from Buzzard in Z.Zhang 2023. lnA is the normalization for the MOR. Different simulatins may 
    give different values. Here we make the assumption that the difference between lnA in an orientatin bin with the mean lnA
    across different simulations yields consistent results. This is because the direction and magnitude of richness 
    boosting effect is consistent with previous literature. 
    */

    //Something wrtong with the contructor. 
    TRIAX_RICHNESS(quad::Interp1D const& lnA)
    :_lnA(lnA)
    {}

    double const lnA_mean = 2.956;
    static std::array<double, 5> constexpr lnA_bin = {2.866, 2.892, 2.916, 2.986, 3.114}; 
    static std::array<double, 5> constexpr mu_bin = {0.1, 0.3, 0.5, 0.7, 0.9};

    explicit TRIAX_RICHNESS(cosmosis::DataBlock& sample)
      : _lnA(mu_bin, lnA_bin)
    {}

    __device__ __host__ double
    operator()(double mu) const
    { //Using the one parameter model for boosting of richness from triaxiality bias using Z.Zhang 2023
      
      //Interpolating lnA as a function of orientatin
      double const result = _lnA.clamp(mu);

      //What we care about is the shifting wrt to the mean of the normalization of the MOR in different orientation bins. 
      double const result2 = result - lnA_mean;
      return result2;
    }
  };
}

#endif
