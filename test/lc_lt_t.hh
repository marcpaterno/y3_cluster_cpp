#ifndef Y3_CLUSTER_LC_LT_T_HH
#define Y3_CLUSTER_LC_LT_T_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "test/primitives.hh"
#include <cmath>
#include <fstream>
#include <iostream>

namespace y3_cluster {
  struct LC_LT_t {
    
    double
    operator()(double lc, double lt, double zt) const
    { 
      /* Define zt_max and l_max large enough to be safe*/
      double constexpr zt_max=2.0;
      double constexpr lt_max=300.0;
      std::array<double, 6> constexpr zt_bins = {0.1,0.15,0.2,0.25,0.3,zt_max};

      std::array<double, 23> constexpr lt_bins = {1.,3.,5.,7.,9.,
                                                12.,15.55555534,20.,24.,26.11111069,
                                                30.,36.66666412,40.,47.22222137,57.77777863,68.33332825,
                                                78.8888855,89.44444275,100.,120.,140.,160.,lt_max};
      
      int constexpr zt_len = zt_bins.size()-1;
      int constexpr lt_len = lt_bins.size()-1;
      int constexpr length = zt_len * lt_len;

      y3_cluster::LC_LT_PARAMS lc_lt_params;

      double tau, mu, sigma, fmsk, fprj;
      tau= mu= sigma= fmsk= fprj=0.0;
      for (std::size_t i = 0; i < zt_len; ++i) {
        for (std::size_t j = 0; j < lt_len; ++j) {
          if ((zt>=zt_bins[i]) && (zt<zt_bins[i+1]) 
                && (lt>=lt_bins[j]) && (lt<lt_bins[j+1])){
            tau = lc_lt_params(i*lt_len+j);
            mu = lc_lt_params(i*lt_len+j+length);
            sigma = lc_lt_params(i*lt_len+j+length*2);
            fmsk = lc_lt_params(i*lt_len+j+length*3);
            fprj = lc_lt_params(i*lt_len+j+length*4);
          }
        }
      }

      return (1.0 - fmsk) * (1.0 - fprj) * y3_cluster::invsqrt2pi() *
               std::exp(-(lc - mu) * (lc - mu) / (2.0 * sigma * sigma)) /
               sigma +
             0.5 * ((1.0 - fmsk) * fprj * tau + fmsk * fprj / lt) *
               std::exp(tau * (2.0 * mu + tau * sigma * sigma - 2.0 * lc) /
                        2.0) *
               std::erfc((mu + tau * sigma * sigma - lc) /
                         (std::sqrt(2.0) * sigma)) +
             fmsk *
               (std::erfc((mu - lc - lt) / (std::sqrt(2.0) * sigma)) -
                std::erfc((mu - lc) / (std::sqrt(2.0) * sigma))) /
               (2.0 * lt) -
             fmsk * fprj *
               (std::exp(-tau * lt) *
                std::exp(tau *
                         (2.0 * mu + tau * sigma * sigma - 2.0 * lc) /
                         2.0) *
                std::erfc((mu + tau * sigma * sigma - lc - lt) /
                          (std::sqrt(2.0) * sigma))) /
               (2.0 * lt);
    }

  };
}

#endif