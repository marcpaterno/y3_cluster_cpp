// Shared weight functors for the radial lensing N-operator drivers.
//
// All lensing profiles are read from the "haloModel" section populated by
// y3_buzzard/haloModelCosmosis.py (see memory/halo_model_cosmosis.md):
//
//   Sigma_nfw (r_sigma, lnM)    -- 1-halo projected surface density
//   Sigma_hh  (r_sigma, z)      -- 2-halo projected surface density (no bias)
//   dSigma_nfw(r_sigma, lnM)    -- 1-halo differential
//   dSigma_hh (r_sigma, z)      -- 2-halo differential (no bias)
//   bias      (lnM,     z)      -- halo bias b(M, z)
//
// Sigma_crit_inv(z_lens) is a source-bin-weighted 1-D function read from
// section "average_sigma_crit_inv" (populated by
// src/modules/average_sigma_crit_inv/average_sigma_crit_inv.py, axes
// "zlense" / "sci_average"). Σ_cr^-1 has no intrinsic R dependence;
// any R dependence in a full weighted pipeline comes from a per-R source
// selection, which the existing smoke module does not model.
//
// Total 1h+2h surface densities (used by the _tot variants):
//   Sigma_tot (R, M, z) = Sigma_nfw(R, M)  + b(M, z) * Sigma_hh (R, z)
//   dSigma_tot(R, M, z) = dSigma_nfw(R, M) + b(M, z) * dSigma_hh(R, z)
//
// Shear uses Sigma_crit_inv to rescale:
//   gamma_t = dSigma * Sigma_crit_inv(z_lens)
//
// The old reduced-shear form g_t = gamma_t / (1 - Sigma * Sigma_crit_inv)
// was retired 2026-05-11: its 1/(1 - x) denominator breaks the additive
// decomposition gamma_t^total = gamma_t^1h + gamma_t^prj that the
// Costanzi-2026 projection observable requires.
//
// Each weight functor exposes:
//   static char const* module_label();
//   explicit Weight(cosmosis::DataBlock&);
//   void set_sample(cosmosis::DataBlock&);
//   double operator()(double R, double lnM, double zt) const;
#ifndef Y3_CLUSTER_CPP_NUM_COUNTS_SEL_LENSING_WEIGHTS_HH
#define Y3_CLUSTER_CPP_NUM_COUNTS_SEL_LENSING_WEIGHTS_HH

#include "cosmosis/datablock/datablock.hh"
#include "utils/interp_1d.hh"
#include "utils/interp_2d.hh"
#include "utils/make_interp_1d.hh"
#include "utils/make_interp_2d.hh"

#include <optional>

namespace y3_cluster_sel_weights {

  // Helper — load Sigma_crit_inv^{-1}(z_lens) 1-D table.
  inline y3_cluster::Interp1D
  load_sigma_crit_inv(cosmosis::DataBlock& s)
  {
    return y3_cluster::make_Interp1D(s, "average_sigma_crit_inv",
                                     "zlense", "sci_average");
  }

  // ---- 1-halo surface density ------------------------------------------------
  struct Sigma1hWeight {
    static char const* module_label() { return "Sigma1hSel"; }
    explicit Sigma1hWeight(cosmosis::DataBlock&) {}
    void
    set_sample(cosmosis::DataBlock& s)
    {
      sigma_nfw_.emplace(
        y3_cluster::make_Interp2D(s, "haloModel", "r_sigma", "lnM", "Sigma_nfw"));
    }
    double
    operator()(double R, double lnM, double /*zt*/) const
    {
      return sigma_nfw_->clamp(R, lnM);
    }
  private:
    std::optional<y3_cluster::Interp2D> sigma_nfw_;
  };


  // ---- Total (1h + b·2h) surface density ------------------------------------
  struct SigmaTotWeight {
    static char const* module_label() { return "SigmaTotSel"; }
    explicit SigmaTotWeight(cosmosis::DataBlock&) {}
    void
    set_sample(cosmosis::DataBlock& s)
    {
      sigma_nfw_.emplace(
        y3_cluster::make_Interp2D(s, "haloModel", "r_sigma", "lnM", "Sigma_nfw"));
      sigma_hh_.emplace(
        y3_cluster::make_Interp2D(s, "haloModel", "r_sigma", "z",   "Sigma_hh"));
      bias_.emplace(
        y3_cluster::make_Interp2D(s, "haloModel", "lnM",     "z",   "bias"));
    }
    double
    operator()(double R, double lnM, double zt) const
    {
      double const b  = bias_->clamp(lnM, zt);
      double const s1 = sigma_nfw_->clamp(R, lnM);
      double const s2 = sigma_hh_ ->clamp(R, zt);
      return s1 + b * s2;
    }
  private:
    std::optional<y3_cluster::Interp2D> sigma_nfw_;
    std::optional<y3_cluster::Interp2D> sigma_hh_;
    std::optional<y3_cluster::Interp2D> bias_;
  };


  // ---- 1-halo differential ΔΣ -----------------------------------------------
  struct DSigma1hWeight {
    static char const* module_label() { return "DSigma1hSel"; }
    explicit DSigma1hWeight(cosmosis::DataBlock&) {}
    void
    set_sample(cosmosis::DataBlock& s)
    {
      dsigma_nfw_.emplace(
        y3_cluster::make_Interp2D(s, "haloModel", "r_sigma", "lnM", "dSigma_nfw"));
    }
    double
    operator()(double R, double lnM, double /*zt*/) const
    {
      return dsigma_nfw_->clamp(R, lnM);
    }
  private:
    std::optional<y3_cluster::Interp2D> dsigma_nfw_;
  };


  // ---- Total (1h + b·2h) differential ΔΣ -----------------------------------
  struct DSigmaTotWeight {
    static char const* module_label() { return "DSigmaTotSel"; }
    explicit DSigmaTotWeight(cosmosis::DataBlock&) {}
    void
    set_sample(cosmosis::DataBlock& s)
    {
      dsigma_nfw_.emplace(
        y3_cluster::make_Interp2D(s, "haloModel", "r_sigma", "lnM", "dSigma_nfw"));
      dsigma_hh_.emplace(
        y3_cluster::make_Interp2D(s, "haloModel", "r_sigma", "z",   "dSigma_hh"));
      bias_.emplace(
        y3_cluster::make_Interp2D(s, "haloModel", "lnM",     "z",   "bias"));
    }
    double
    operator()(double R, double lnM, double zt) const
    {
      double const b  = bias_->clamp(lnM, zt);
      double const d1 = dsigma_nfw_->clamp(R, lnM);
      double const d2 = dsigma_hh_ ->clamp(R, zt);
      return d1 + b * d2;
    }
  private:
    std::optional<y3_cluster::Interp2D> dsigma_nfw_;
    std::optional<y3_cluster::Interp2D> dsigma_hh_;
    std::optional<y3_cluster::Interp2D> bias_;
  };


  // ---- 1-halo shear γ_t = ΔΣ_nfw · Σ_crit^-1(z_lens) -----------------------
  struct Shear1hWeight {
    static char const* module_label() { return "Shear1hSel"; }
    explicit Shear1hWeight(cosmosis::DataBlock&) {}
    void
    set_sample(cosmosis::DataBlock& s)
    {
      dsigma_nfw_.emplace(
        y3_cluster::make_Interp2D(s, "haloModel", "r_sigma", "lnM", "dSigma_nfw"));
      sigma_crit_inv_.emplace(load_sigma_crit_inv(s));
    }
    double
    operator()(double R, double lnM, double zt) const
    {
      double const d = dsigma_nfw_->clamp(R, lnM);
      double const sci = sigma_crit_inv_->clamp(zt);
      return d * sci;
    }
  private:
    std::optional<y3_cluster::Interp2D> dsigma_nfw_;
    std::optional<y3_cluster::Interp1D> sigma_crit_inv_;
  };


  // ---- Total shear γ_t = ΔΣ_tot · Σ_crit^-1(z_lens) ------------------------
  struct ShearTotWeight {
    static char const* module_label() { return "ShearTotSel"; }
    explicit ShearTotWeight(cosmosis::DataBlock&) {}
    void
    set_sample(cosmosis::DataBlock& s)
    {
      dsigma_nfw_.emplace(
        y3_cluster::make_Interp2D(s, "haloModel", "r_sigma", "lnM", "dSigma_nfw"));
      dsigma_hh_.emplace(
        y3_cluster::make_Interp2D(s, "haloModel", "r_sigma", "z",   "dSigma_hh"));
      bias_.emplace(
        y3_cluster::make_Interp2D(s, "haloModel", "lnM",     "z",   "bias"));
      sigma_crit_inv_.emplace(load_sigma_crit_inv(s));
    }
    double
    operator()(double R, double lnM, double zt) const
    {
      double const b = bias_->clamp(lnM, zt);
      double const d = dsigma_nfw_->clamp(R, lnM) +
                       b * dsigma_hh_->clamp(R, zt);
      double const sci = sigma_crit_inv_->clamp(zt);
      return d * sci;
    }
  private:
    std::optional<y3_cluster::Interp2D> dsigma_nfw_;
    std::optional<y3_cluster::Interp2D> dsigma_hh_;
    std::optional<y3_cluster::Interp2D> bias_;
    std::optional<y3_cluster::Interp1D> sigma_crit_inv_;
  };


  // ReducedShear{1h,Tot}Weight structs retired 2026-05-11.  We now use
  // the linear γ_t = ΔΣ · Σ_crit^-1 form (Shear{1h,Tot}Weight above);
  // the 1h + prj observable is assembled additively in the likelihood,
  // which was not valid under the old 1/(1 - ΣΣ_crit^-1) denominator.

} // namespace y3_cluster_sel_weights

#endif
