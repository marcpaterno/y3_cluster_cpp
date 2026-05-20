// RmSelFunction_t -- closed-form bin-integrated richness kernel K_i(ltr, z)
// for the DES Y3 selection function.
//
// Implements Eq. (34) (Ki_final_expanded) of
// $PSCRATCH/github/RichnessSelection/docs/richness_selection_function.tex,
// in compact form Eq. (33) (Ki_final):
//
//   K_i(ltr, z) = (1 - fprj) [Phi((lmax - mu) / sigma)
//                             - Phi((lmin - mu) / sigma)]
//               + fprj       [F_EMG(lmax) - F_EMG(lmin)]
//
// with mu = ltr + Delta_mu and the EMG parameters (Delta_mu, sigma, tau,
// fprj) read as functions of (ltr, z) from the Costanzi
// `prj_params_DESY3_lss_lin_dep_getdist_v1.txt` calibration via
// PlobLtrEMG_t (datablock section "plob_ltr_params").
//
// The numerics live in RichnessKernel_t (richness_kernel_t.hh): this class
// is a thin (bin_index -> bin edges) + (datablock -> PlobLtrEMG_t)
// adapter so the Costanzi-2026 integrands can address the kernel by bin.
//
// Mirrors richness_selection/selection_function/kernels.py::K_i.
#ifndef Y3_CLUSTER_PROJECTION_Y3_B_I_T_HH
#define Y3_CLUSTER_PROJECTION_Y3_B_I_T_HH

#include "cosmosis/datablock/datablock.hh"
#include "models/plob_ltr_emg_t.hh"
#include "models/richness_kernel_t.hh"

#include <array>
#include <limits>

namespace y3_cluster {

  class RmSelFunction_t {
  public:
    // DES Y3 richness bin edges [20, 30, 45, 60, inf).
    static constexpr std::array<double, 5> bin_edges{
      {20.0, 30.0, 45.0, 60.0, std::numeric_limits<double>::infinity()}};
    static constexpr int n_bins = 4;

    RmSelFunction_t(cosmosis::DataBlock& sample, int bin_index)
      : _plob(sample)
      , _kernel(bin_edges[bin_index], bin_edges[bin_index + 1])
    {}

    double
    operator()(double lt, double z) const
    {
      return _kernel(lt, z, _plob);
    }

  private:
    PlobLtrEMG_t      _plob;
    RichnessKernel_t  _kernel;
  };

} // namespace y3_cluster

#endif
