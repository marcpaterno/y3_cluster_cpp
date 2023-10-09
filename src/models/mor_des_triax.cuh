#ifndef Y3_CLUSTER_CPP_MOR_DES_TRIAX_CUH
#define Y3_CLUSTER_CPP_MOR_DES_TRIAX_CUH

#include "cosmosis/datablock/datablock.hh"
#include "common/cuda/Interp2D.cuh"
#include "utils/primitives.cuh"
#include "models/triax_richness.cuh"

#include <array>
#include <cmath>

/*
Z.Zhang. From Z.Zhang+22 paper for the triaxiaity model. 
*/

namespace y3_cuda {

  class MOR_DES_TRIAX {
    double _lnA;
    double _B;
    double _C;
    double _sigma_intr;
    double _m_pivot;
    double _z_pivot;

  private:
    std::optional<y3_cuda::TRIAX_RICHNESS> dlnA_cosi; 

  public:
    MOR_DES_TRIAX(double lnA,
              double B,
              double C,
              double sigma_intr,
              double m_pivot,
              double z_pivot)
    {}

    explicit MOR_DES_TRIAX(cosmosis::DataBlock& sample):
      _lnA(sample.view<double>("mor_loglinear", "lnA")),
      _B(sample.view<double>("mor_loglinear", "B")), 
      _C(sample.view<double>("mor_loglinear", "C")), 
      _sigma_intr(sample.view<double>("mor_loglinear", "sigma_intr")),
      _m_pivot(sample.view<double>("mor_loglinear", "m_pivot")),
      _z_pivot(sample.view<double>("mor_loglinear", "mor_pivot"))
    {}

    __device__ __host__ double
    operator()(double lt, double lnM, double zt, double mu) const
    {
      //Mean MOR
      double const lnl = _lnA+ _B * (lnM - std::log(_m_pivot)) + _C*std::log((1+zt)/(1+_z_pivot)); 
      
      //Boosting from triaxiality. TODO: think about polynomial solution. 
      double const delta_lnA_cosi = (*dlnA_cosi)(mu);
      double const lnl_tot = lnl + delta_lnA_cosi;
      
      double const x = lt - lnl_tot;
      double const _sigma = std::sqrt(std::pow(_sigma_intr,2) + std::exp(lnl_tot)-1/std::exp(2*lnl_tot)); 

      //!!TODO!! Think about finding analytic integration solution as a function of mu.
      return y3_cuda::gaussian(x, 0.0, _sigma);
    }
  };
}

#endif
