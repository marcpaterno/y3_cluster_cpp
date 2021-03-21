#ifndef Y3_CLUSTER_CPP_MOR_DES_T_HH
#define Y3_CLUSTER_CPP_MOR_DES_T_HH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "ln_mez_power_law.hh"
#include "utils/interp_2d.hh"
#include "utils/make_interp_2d.hh"
#include "utils/mz_power_law.hh"

#include <cmath>

namespace y3_cluster {

  using doubles = std::vector<double>;

  // This MOR is a bivariate gaussian in lamba^true and ln(zeta) given
  // ln(M200m) and z^true.
  //      <lam^true> = 1 + <lam^sat>
  //      <lam^sat> = ((M200m - Mmin) / (M1 - Mmin))^alpha * ((1+z) /
  //      (1+zp))^epsilon Variance(lam^true) = <lam^sat> + (<lam^sat> *
  //      sigma_intr)^2 <ln zeta> = ln Asz + Bsz * (lnM500c - lnMpsz) + Csz *
  //      (lnE(z) - lnE(zpsz)) Variance(ln zeta) = sigma_lnzeta^2
  class MOR_Y3xSPT_t {
  public:
    MOR_Y3xSPT_t(double Mmin,
                 double M1,
                 double alpha,
                 double epsilon,
                 double sigintr,
                 double Asz,
                 double Bsz,
                 double Csz,
                 double siglnzeta,
                 double r,
                 double zplam,
                 double lnMpsz,
                 double zpsz,
                 double omega_m,
                 double omega_l,
                 double omega_k,
                 Interp2D const& mass_conversion)
      : Mmin_(Mmin)
      , M1_(M1)
      , alpha_(alpha)
      , epsilon_(epsilon)
      , sigintr_(sigintr)
      , Asz_(Asz)
      , Bsz_(Bsz)
      , Csz_(Csz)
      , siglnzeta_(siglnzeta)
      , r_(r)
      , zplam_(zplam)
      , lnMpsz_(lnMpsz)
      , zpsz_(zpsz)
      , omega_m_(omega_m)
      , omega_l_(omega_l)
      , omega_k_(omega_k)
      , zetamrel_(Asz_,
                  Bsz_,
                  Csz_,
                  lnMpsz_,
                  zpsz_,
                  omega_m_,
                  omega_l_,
                  omega_k_)
      , mass_conversion_(mass_conversion)
    {}

    MOR_Y3xSPT_t(cosmosis::DataBlock& sample)
      : Mmin_(sample.view<double>("abundance_integral", "Mmin"))
      , M1_(sample.view<double>("abundance_integral", "M1"))
      , alpha_(sample.view<double>("abundance_integral", "alpha"))
      , epsilon_(sample.view<double>("abundance_integral", "epsilon"))
      , sigintr_(sample.view<double>("abundance_integral", "sig0lam"))
      , Asz_(sample.view<double>("abundance_integral", "Asz"))
      , Bsz_(sample.view<double>("abundance_integral", "Bsz"))
      , Csz_(sample.view<double>("abundance_integral", "Csz"))
      , siglnzeta_(sample.view<double>("abundance_integral", "sig0zeta"))
      , r_(sample.view<double>("abundance_integral", "r"))
      , zplam_(sample.view<double>("abundance_integral", "zplam"))
      , lnMpsz_(sample.view<double>("abundance_integral", "lnMpsz"))
      , zpsz_(sample.view<double>("abundance_integral", "zpsz"))
      , omega_m_(sample.view<double>("cosmological_parameters", "omega_m"))
      , omega_l_(get_datablock<double>(sample,
                                       "cosmological_parameters",
                                       "omega_lambda"))
      , omega_k_(sample.view<double>("cosmological_parameters", "omega_k"))
      , zetamrel_(Asz_,
                  Bsz_,
                  Csz_,
                  lnMpsz_,
                  zpsz_,
                  omega_m_,
                  omega_l_,
                  omega_k_)
      , mass_conversion_(
          make_Interp2D(sample, "mass_conversion", "lnm200m", "z", "lnm500c"))
    {}

    double
    operator()(double lamtrue,
               double zeta,
               double ztrue,
               double lnM200m,
               double gamma_field) const
    {
      // Richness scaling calculations
      // Can precalc some of this in constructor...
      // Or make it into a callable model for independent testing
      double const mean_lamsat =
        pow((std::exp(lnM200m) - Mmin_) / (M1_ - Mmin_), alpha_) *
        pow((1.0 + ztrue) / (1.0 + zplam_), epsilon_);
      double const mean_lamtrue = 1.0 + mean_lamsat;
      double const var_lam =
        mean_lamsat + mean_lamsat * mean_lamsat * sigintr_ * sigintr_;

      // ln zeta scaling calculations
      double const lnM500c = mass_conversion_(lnM200m, ztrue);
      double const mean_lnzeta = zetamrel_(lnM500c, ztrue, gamma_field);
      double const var_lnzeta = siglnzeta_ * siglnzeta_;

      // Compute repeated quantities
      double const delta_lamtrue = lamtrue - mean_lamtrue;
      double const delta_lnzeta = std::log(zeta) - mean_lnzeta;
      double const one_over_det_cov =
        1.0 / (var_lnzeta * var_lam * (1.0 - r_ * r_));

      double prefactor = 0.15915494309 * std::sqrt(one_over_det_cov) / zeta;
      double term1 = delta_lamtrue * delta_lamtrue * var_lnzeta;
      double term2 = -2.0 * r_ * std::sqrt(var_lam * var_lnzeta) *
                     delta_lamtrue * delta_lnzeta;
      double term3 = delta_lnzeta * delta_lnzeta * var_lam;

      return prefactor *
             std::exp(-0.5 * one_over_det_cov * (term1 + term2 + term3));
    }

  private:
    double const Mmin_, M1_, alpha_, epsilon_, sigintr_;
    double const Asz_, Bsz_, Csz_, siglnzeta_, r_;
    double const zplam_, lnMpsz_, zpsz_;
    double const omega_m_, omega_l_, omega_k_;
    y3_cluster::ln_mez_power_law zetamrel_;
    Interp2D mass_conversion_;
  };
}

#endif
