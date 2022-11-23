#ifndef Y3_CUDA_LC_LT_T_CUH
#define Y3_CUDA_LC_LT_T_CUH

// jta conversion project to cuda Nov 2022

#include "cosmosis/datablock/datablock.hh"
#include "cuda/pagani/quad/GPUquad/Interp2D.cuh"
#include "utils/primitives.cuh"

// Implementing Costanzi, Rozo, Rykoff et al 2018 equation 15
// which is a P(\lambda^obs | \lambda^true, z) form
// incorporating background, projection, and percolation
// 
// In Yuanyaun's section 9 summary of the Lighthouse document
// https://www.overleaf.com/project/5c378b07f882d02f5b8c90e2
// Costanzi's eq 15 is labeled P(\lambda_c|\lambda_t) (eq 51).

namespace y3_cuda {
  struct LC_LT_t {

    // CUDA/C++ does not support static data members of classes or structs
    quad::Interp2D tau_interp;
    quad::Interp2D mu_interp;
    quad::Interp2D sigma_interp;
    quad::Interp2D fmsk_interp;
    quad::Interp2D fprj_interp;

    // The default constructor is implemented in the .cu file
    LC_LT_t() ;

    // The constructor from a datablock is required, but does the same thing
    // as the default constructor. So we delegate the work to the default
    // constructor.
    explicit LC_LT_t(cosmosis::DataBlock const&) : LC_LT_t() { }

    __device__ __host__ double
    operator()(double lc, double lt, double zt) const
    {
      const auto tau = tau_interp.clamp(lt, zt);
      const auto mu = mu_interp.clamp(lt, zt);
      const auto sigma = sigma_interp.clamp(lt, zt);
      const auto fmsk = fmsk_interp.clamp(lt, zt);
      const auto fprj = fprj_interp.clamp(lt, zt);

      const auto exptau =
        std::exp(tau * (2.0 * mu + tau * sigma * sigma - 2.0 * lc) / 2.0);
      const auto root_two_sigma = std::sqrt(2.0) * sigma;
      const auto mu_tau_sig_sqr = mu + tau * sigma * sigma;

      // Helper function for common pattern
      const auto erfc_scaled = [root_two_sigma](double a, double b) {
        return std::erfc((a - b) / root_two_sigma);
      };

      // eq. (33)     https://www.overleaf.com/project/5c378b07f882d02f5b8c90e2
      // But also eq 51 in Yuanyuan's clean up section in that paper
      // and
      // Implementing Costanzi, Rozo, Rykoff et al 2019 equation 15
      // which is a P(\lambda^obs | \lambda^true, z) form
      return (1.0 - fmsk) * (1.0 - fprj) * y3_cuda::gaussian(lc, mu, sigma) +
             0.5 * ((1.0 - fmsk) * fprj * tau + fmsk * fprj / lt) * exptau *
               erfc_scaled(mu_tau_sig_sqr, lc) +
             0.5 * fmsk / lt *
               (erfc_scaled(lc, mu) - erfc_scaled(lc + lt, mu)) -
             0.5 * fmsk * fprj / lt *
               (std::exp(-tau * lt) * exptau *
                erfc_scaled(mu_tau_sig_sqr, lc + lt));
    }
  };
}

#endif
