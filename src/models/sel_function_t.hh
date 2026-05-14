// SelFunction_t — per-bin Interp2D reader for the richness-selection
// function S_ij(lnM, z) precomputed by sel_function.py.
//
// Data layout (single packed tensor, shared axes):
//   sel_function/lnM      (n_lnm,)
//   sel_function/z        (n_z,)
//   sel_function/S_stack  (n_bins, n_z, n_lnm)
//
// The feeder writes all 12 bin tables as one contiguous ndarray for
// speed (one put_val instead of 12). At construction we slice the
// requested bin's (n_z, n_lnm) plane and pass it to Interp2D.
#ifndef Y3_CLUSTER_CPP_SEL_FUNCTION_T_HH
#define Y3_CLUSTER_CPP_SEL_FUNCTION_T_HH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "utils/interp_2d.hh"

#include <algorithm>
#include <string>
#include <vector>

namespace y3_cluster {

  class SelFunction_t {
  public:
    using doubles = std::vector<double>;

    SelFunction_t(cosmosis::DataBlock& sample, int bin_index)
      : _p(_make_interp(sample, bin_index))
    {
      auto const& lnm = sample.view<doubles>("sel_function", "lnM");
      auto const& z   = sample.view<doubles>("sel_function", "z");
      auto [lnm_min, lnm_max] = std::minmax_element(lnm.begin(), lnm.end());
      auto [z_min,   z_max]   = std::minmax_element(z.begin(),   z.end());
      _lnm_lo = *lnm_min; _lnm_hi = *lnm_max;
      _z_lo   = *z_min;   _z_hi   = *z_max;
    }

    double
    operator()(double lnM, double zt) const
    {
      if (lnM < _lnm_lo || lnM > _lnm_hi || zt < _z_lo || zt > _z_hi)
        return 0.0;
      return _p.clamp(lnM, zt);
    }

  private:
    Interp2D _p;
    double _lnm_lo{};
    double _lnm_hi{};
    double _z_lo{};
    double _z_hi{};

    static Interp2D
    _make_interp(cosmosis::DataBlock& sample, int bin_index)
    {
      auto const& lnm = sample.view<doubles>("sel_function", "lnM");
      auto const& z   = sample.view<doubles>("sel_function", "z");
      auto const& nd  = sample.view<cosmosis::ndarray<double>>(
        "sel_function", "S_stack");
      auto ext = nd.extents();
      if (ext.size() != 3) {
        throw std::runtime_error(
          "sel_function/S_stack must be 3-D (n_bins, n_z, n_lnm)");
      }
      std::size_t const nb  = ext[0];
      std::size_t const nz  = ext[1];
      std::size_t const nmx = ext[2];
      if (bin_index < 0 || static_cast<std::size_t>(bin_index) >= nb) {
        throw std::runtime_error("sel_function: bin_index out of range");
      }
      if (nz != z.size() || nmx != lnm.size()) {
        throw std::runtime_error("sel_function: S_stack extents mismatch");
      }
      std::size_t const plane = nz * nmx;
      std::size_t const offset = static_cast<std::size_t>(bin_index) * plane;
      // Copy the requested plane into a std::vector and hand it to Interp2D.
      std::vector<double> vals(nd.begin() + offset,
                                nd.begin() + offset + plane);
      return Interp2D(lnm, z, vals);
    }
  };

} // namespace y3_cluster

#endif
