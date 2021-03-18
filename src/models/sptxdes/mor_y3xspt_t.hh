#ifndef Y3_CLUSTER_CPP_MOR_DES_T_HH
#define Y3_CLUSTER_CPP_MOR_DES_T_HH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "ln_mez_power_law.hh"
#include "utils/datablock_reader.hh"
#include "utils/interp_2d.hh"
#include "utils/mz_power_law.hh"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

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
                 std::shared_ptr<Interp2D const> mass_conversion)
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
      : Mmin_(get_datablock<double>(sample, "abundance_integral", "Mmin"))
      , M1_(get_datablock<double>(sample, "abundance_integral", "M1"))
      , alpha_(get_datablock<double>(sample, "abundance_integral", "alpha"))
      , epsilon_(get_datablock<double>(sample, "abundance_integral", "epsilon"))
      , sigintr_(get_datablock<double>(sample, "abundance_integral", "sig0lam"))
      , Asz_(get_datablock<double>(sample, "abundance_integral", "Asz"))
      , Bsz_(get_datablock<double>(sample, "abundance_integral", "Bsz"))
      , Csz_(get_datablock<double>(sample, "abundance_integral", "Csz"))
      , siglnzeta_(
          get_datablock<double>(sample, "abundance_integral", "sig0zeta"))
      , r_(get_datablock<double>(sample, "abundance_integral", "r"))
      , zplam_(get_datablock<double>(sample, "abundance_integral", "zplam"))
      , lnMpsz_(get_datablock<double>(sample, "abundance_integral", "lnMpsz"))
      , zpsz_(get_datablock<double>(sample, "abundance_integral", "zpsz"))
      , omega_m_(
          get_datablock<double>(sample, "cosmological_parameters", "omega_m"))
      , omega_l_(get_datablock<double>(sample,
                                       "cosmological_parameters",
                                       "omega_lambda"))
      , omega_k_(
          get_datablock<double>(sample, "cosmological_parameters", "omega_k"))
      , zetamrel_(Asz_,
                  Bsz_,
                  Csz_,
                  lnMpsz_,
                  zpsz_,
                  omega_m_,
                  omega_l_,
                  omega_k_)
      , mass_conversion_(std::make_shared<Interp2D const>(
          get_datablock<doubles>(sample, "mass_conversion", "lnm200m"),
          get_datablock<doubles>(sample, "mass_conversion", "z"),
          get_datablock<cosmosis::ndarray<double>>(sample,
                                                   "mass_conversion",
                                                   "lnm500c")))
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
      double const lnM500c = mass_conversion_->eval(lnM200m, ztrue);
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
    std::shared_ptr<Interp2D const> mass_conversion_;
  };
}

#endif
