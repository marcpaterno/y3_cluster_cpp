// K_i(ltr, z) and K_j(ztr) — closed-form bin-integrated richness and
// redshift kernels for the DES Y3 selection function.
//
// Derivation in docs/richness_selection_function.tex §3 (K_i) + §1 (K_j)
// of the RichnessSelection repo. Matches
// src/richness_selection/selection_function/kernels.py line for line.
//
//   K_i(ltr, z) = (1 - fprj) [Phi(z_max) - Phi(z_min)]
//               + fprj       [F_EMG(lmax) - F_EMG(lmin)]
//
// with z_{min,max} = (lam_{min,max} - mu) / sigma and F_EMG the EMG CDF
// (eq. 19 of the doc), written through scipy's erfcx(t) = exp(t^2) erfc(t)
// so both tails stay finite in double precision.
//
// K_i is parameterised by (lob_min, lob_max) bin edges and an EMG-parameter
// provider `PlobT` (template) that exposes mu/sigma/tau/fprj(ltr, z).
// Typically PlobT = y3_cluster::PlobLtrEMG_t; templated so tests can swap in
// a stub.
#ifndef Y3_CLUSTER_CPP_RICHNESS_KERNEL_T_HH
#define Y3_CLUSTER_CPP_RICHNESS_KERNEL_T_HH

#include <algorithm>
#include <cmath>

namespace y3_cluster {

  namespace rk_detail {
    inline constexpr double SQRT2 = 1.41421356237309504880;

    inline double phi(double x) {
      return 0.5 * (1.0 + std::erf(x / SQRT2));
    }

    // erfcx(t) = exp(t^2) * erfc(t). Stable asymptotic expansion for
    // t >= ERFCX_CUTOFF; uses exp(t^2)*erfc(t) directly otherwise.
    // Matches scipy.special.erfcx to ~1e-15 at t < 30.
    inline double erfcx(double t) {
      constexpr double ERFCX_CUTOFF = 25.0;
      if (t < ERFCX_CUTOFF) {
        return std::exp(t * t) * std::erfc(t);
      }
      // Asymptotic series 1/(sqrt(pi)*t) * (1 - 1/(2t^2) + 3/(4t^4) - ...)
      double const inv = 1.0 / t;
      double const inv2 = inv * inv;
      return inv / std::sqrt(M_PI) *
             (1.0 - 0.5 * inv2 * (1.0 - 1.5 * inv2 * (1.0 - 2.5 * inv2)));
    }

    // F_EMG(x; mu, sigma, tau) via erfcx — matches kernels.py::_F_EMG_vec.
    //
    //   F_EMG(x) = Phi(z) - tail(x, mu, sigma, tau)
    //   with z = (x - mu) / sigma, u = (tau*sigma - z) / sqrt(2)
    //   tail = 0.5 * erfcx(|u|) * exp(-0.5 z^2)        (u >= 0)
    //        = exp(A) - 0.5 * erfcx(|u|) * exp(-0.5 z^2)  (u <  0)
    //   A = -tau (x - mu) + 0.5 tau^2 sigma^2
    inline double
    F_EMG(double x, double mu, double sigma, double tau)
    {
      double const z = (x - mu) / sigma;
      double const u = (tau * sigma - z) / SQRT2;
      bool const neg = u < 0.0;
      double const abs_u = neg ? -u : u;
      double const exp_mz2 = std::exp(-0.5 * z * z);
      double const tail_base = 0.5 * erfcx(abs_u) * exp_mz2;
      double tail;
      if (neg) {
        double const A = -tau * (x - mu) + 0.5 * tau * tau * sigma * sigma;
        tail = std::exp(A) - tail_base;
      } else {
        tail = tail_base;
      }
      return std::clamp(phi(z) - tail, 0.0, 1.0);
    }

  } // namespace rk_detail


  // Bin-integrated richness kernel K_i(ltr, z). Constructed with bin edges;
  // evaluated against the EMG-parameter provider at each (ltr, z).
  class RichnessKernel_t {
  public:
    RichnessKernel_t(double lob_min, double lob_max)
      : _lob_min(lob_min), _lob_max(lob_max)
    {}

    // Evaluate K_i at a single (ltr, z), using plob's scalar mu/sigma/tau/fprj.
    template <typename PlobT>
    double
    operator()(double ltr, double z, PlobT const& plob) const
    {
      double const mu    = plob.mu   (ltr, z);
      double const sigma = plob.sigma(ltr, z);
      double const tau   = plob.tau  (ltr, z);
      double const fprj  = std::min(1.0, plob.fprj (ltr, z));

      double const gauss_piece =
        rk_detail::phi((_lob_max - mu) / sigma) -
        rk_detail::phi((_lob_min - mu) / sigma);

      double const emg_piece =
        rk_detail::F_EMG(_lob_max, mu, sigma, tau) -
        rk_detail::F_EMG(_lob_min, mu, sigma, tau);

      return (1.0 - fprj) * gauss_piece + fprj * emg_piece;
    }

    double lob_min() const { return _lob_min; }
    double lob_max() const { return _lob_max; }

  private:
    double _lob_min;
    double _lob_max;
  };


  // K_j(ztr) = Phi((zob_max - ztr) / sigma_z) - Phi((zob_min - ztr) / sigma_z)
  //
  // Gaussian photo-z kernel — sigma_z is the per-richness-bin photo-z
  // scatter scalar (doc Eq. 12).
  inline double
  richness_zkernel(double ztr, double zob_min, double zob_max, double sigma_z)
  {
    return rk_detail::phi((zob_max - ztr) / sigma_z) -
           rk_detail::phi((zob_min - ztr) / sigma_z);
  }

} // namespace y3_cluster

#endif
