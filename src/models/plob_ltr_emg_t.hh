// P(lob | ltr, z) EMG kernel parameters — DES Y3 lss_lin_dep calibration.
//
// Mirrors src/richness_selection/plob_ltr.py in the RichnessSelection Python
// package (data/prj_params_DESY3_lss_lin_dep_getdist_v1.txt). The four EMG
// parameters (mu, sigma, tau, fprj) are each closed-form functions of ltr
// whose coefficients (a(z), b(z)) are Interp1D splines over 15 z-nodes in
// [0.10, 0.80]:
//
//   mu(ltr, z)    = a_mu(z)   + b_mu(z)   * ltr
//   sigma(ltr, z) = b_sig(z)  * ltr ** a_sig(z)
//   tau(ltr, z)   = b_tau(z)  / ltr ** a_tau(z)
//   fprj(ltr, z)  = b_fprj(z) / (1 + exp(-ltr)) ** a_fprj(z)
//
// These are the EMG parameters for the Costanzi P(lob|ltr,z) kernel
// (projection effects + Gaussian noise; see kernels.py::F_EMG and
// docs/richness_selection_function.tex §3).
//
// Datablock section "plob_ltr_params" is expected to carry:
//   z, a_tau, b_tau, a_mu, b_mu, a_sig, b_sig, a_fprj, b_fprj
// (length-15 double arrays). The companion Python module
// y3_buzzard/prj_params.py (PrjParams) publishes these once at
// pipeline setup; tests populate them inline.
#ifndef Y3_CLUSTER_CPP_PLOB_LTR_EMG_T_HH
#define Y3_CLUSTER_CPP_PLOB_LTR_EMG_T_HH

#include "cosmosis/datablock/datablock.hh"
#include "utils/interp_1d.hh"
#include "utils/make_interp_1d.hh"

#include <algorithm>
#include <cassert>
#include <cmath>

namespace y3_cluster {

  class PlobLtrEMG_t {
  public:
    explicit PlobLtrEMG_t(cosmosis::DataBlock& sample)
      : _a_tau (make_Interp1D(sample, "plob_ltr_params", "z", "a_tau"))
      , _b_tau (make_Interp1D(sample, "plob_ltr_params", "z", "b_tau"))
      , _a_mu  (make_Interp1D(sample, "plob_ltr_params", "z", "a_mu"))
      , _b_mu  (make_Interp1D(sample, "plob_ltr_params", "z", "b_mu"))
      , _a_sig (make_Interp1D(sample, "plob_ltr_params", "z", "a_sig"))
      , _b_sig (make_Interp1D(sample, "plob_ltr_params", "z", "b_sig"))
      , _a_fprj(make_Interp1D(sample, "plob_ltr_params", "z", "a_fprj"))
      , _b_fprj(make_Interp1D(sample, "plob_ltr_params", "z", "b_fprj"))
    {}

    double mu   (double ltr, double z) const {
      return _a_mu.clamp(z) + _b_mu.clamp(z) * ltr;
    }
    double sigma(double ltr, double z) const {
      assert(ltr > 0.0);
      return _b_sig.clamp(z) * std::pow(ltr, _a_sig.clamp(z));
    }
    double tau  (double ltr, double z) const {
      assert(ltr > 0.0);
      return _b_tau.clamp(z) / std::pow(ltr, _a_tau.clamp(z));
    }
    double fprj (double ltr, double z) const {
      double const a = _a_fprj.clamp(z);
      double const b = _b_fprj.clamp(z);
      double const denom = std::pow(1.0 + std::exp(-ltr), a);
      return std::min(1.0, b / denom);
    }

    // Full P(lob | ltr, z) EMG kernel (for debugging / tests); production
    // code does the bin integral analytically via RichnessKernel_t.
    double
    operator()(double lob, double ltr, double z) const
    {
      double const m = mu(ltr, z);
      double const s = sigma(ltr, z);
      double const t = tau(ltr, z);
      double const f = fprj(ltr, z);

      double const u = (lob - m) / s;
      double const gauss = std::exp(-0.5 * u * u) /
                           (s * std::sqrt(2.0 * M_PI));
      double const exp_arg  = 0.5 * t * (2.0 * m + t * s * s - 2.0 * lob);
      double const erfc_arg = (m + t * s * s - lob) / (std::sqrt(2.0) * s);
      double const emg = 0.5 * t * std::exp(exp_arg) * std::erfc(erfc_arg);

      return (1.0 - f) * gauss + f * emg;
    }

  private:
    Interp1D _a_tau,  _b_tau;
    Interp1D _a_mu,   _b_mu;
    Interp1D _a_sig,  _b_sig;
    Interp1D _a_fprj, _b_fprj;
  };

} // namespace y3_cluster

#endif
