// B_i(λ^true, z) — closed-form bin integral of the DES Y3 projection kernel
//
//   P(λ^ob | λ^true, z) = (1 - f_prj) · δ(λ^ob - λ^true)
//                       + f_prj · τ · exp(-τ (λ^ob - λ^true)) · Θ(λ^ob ≥ λ^true)
//
// where (f_prj, τ) = (f_prj, τ)(λ^true, z) are the frozen Y3 tables.
//
//   B_i(lt, z) = (1 - f_prj) · 1[L_i ≤ lt < U_i]
//              + f_prj · [exp(-τ · max(0, L_i - lt))
//                         - exp(-τ · max(0, U_i - lt))]
//
// For the last bin (U_i = ∞) the second exponential vanishes.
//
// Matches `kernels/projection.py:111-154` in the LambdaSelectionEmulator repo.
#ifndef Y3_CLUSTER_PROJECTION_Y3_B_I_T_HH
#define Y3_CLUSTER_PROJECTION_Y3_B_I_T_HH

#include "utils/interp_2d.hh"
#include "models/lc_lt_projection_y3.hh"

#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace y3_cluster {

  class PROJ_Y3_B_i_t {
  public:
    // Bin edges [20, 30, 45, 60, ∞).
    static constexpr std::array<double, 5> bin_edges{
      {20.0, 30.0, 45.0, 60.0, std::numeric_limits<double>::infinity()}};
    static constexpr int n_bins = 4;

    explicit PROJ_Y3_B_i_t(int bin_index)
      : _fprj(_lt_vec(), _zt_vec(), _fprj_vec())
      , _tau(_lt_vec(), _zt_vec(), _tau_vec())
      , _Llo(bin_edges[bin_index])
      , _Lhi(bin_edges[bin_index + 1])
      , _upper_is_inf(std::isinf(bin_edges[bin_index + 1]))
    {}

    double
    operator()(double lt, double zt) const
    {
      double const fp = std::min(std::max(_fprj.clamp(lt, zt), 0.0), 1.0);
      double const tt = std::max(_tau.clamp(lt, zt), 1e-8);

      double const indicator = (lt >= _Llo && lt < _Lhi) ? 1.0 : 0.0;

      double const dL = std::max(0.0, _Llo - lt);
      double const expL = std::exp(-tt * dL);

      double expU = 0.0;
      if (!_upper_is_inf) {
        double const dU = std::max(0.0, _Lhi - lt);
        expU = std::exp(-tt * dU);
      }

      return (1.0 - fp) * indicator + fp * (expL - expU);
    }

  private:
    Interp2D _fprj;
    Interp2D _tau;
    double _Llo;
    double _Lhi;
    bool _upper_is_inf;

    static std::vector<double>
    _lt_vec()
    {
      return std::vector<double>(projection_y3::lt_bins.begin(),
                                 projection_y3::lt_bins.end());
    }
    static std::vector<double>
    _zt_vec()
    {
      return std::vector<double>(projection_y3::zt_bins.begin(),
                                 projection_y3::zt_bins.end());
    }
    static std::vector<double>
    _fprj_vec()
    {
      return std::vector<double>(projection_y3::fprj_arr.begin(),
                                 projection_y3::fprj_arr.end());
    }
    static std::vector<double>
    _tau_vec()
    {
      return std::vector<double>(projection_y3::tau_arr.begin(),
                                 projection_y3::tau_arr.end());
    }
  };

} // namespace y3_cluster

#endif
