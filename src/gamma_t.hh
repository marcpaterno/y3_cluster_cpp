#ifndef Y3_CLUSTER_GAMMA_T_HH
#define Y3_CLUSTER_GAMMA_T_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include <integration_range.hh>
#include <transform.hh>

#include <algorithm>
#include <array>
#include <cubacpp/cubacpp.hh>
#include <cmath>
#include <iostream>

// This class template is based
// on https://www.overleaf.com/13697016cyvvqqfchfbg#/52989522/, and the example
// provided by Spencer Everett.
namespace y3_cluster {

// Forward declaration so we can declare it friend
template <std::size_t NRADII>
struct Gamma_T_Integrated_Bin_Result;

template <typename MODELS, std::size_t NRADII, std::size_t NRICHNESS = 1, std::size_t NREDSHIFT = 1>
class Gamma_T_Integrand {
  friend struct Gamma_T_Integrated_Bin_Result<NRADII>;
private:
  double fcen_;
  double msci_;

  typename MODELS::MOR mor;
  typename MODELS::LO_LC lo_lc;
  typename MODELS::LC_LT lc_lt;
  typename MODELS::ZO_ZT zo_zt;
  typename MODELS::ROFFSET roffset;
  typename MODELS::T_CEN T_cen;
  typename MODELS::T_MIS T_mis;
  typename MODELS::A_CEN A_cen;
  typename MODELS::A_MIS A_mis;
  typename MODELS::HMF hmf;
  typename MODELS::DEL_SIG del_sig;
  typename MODELS::DV_DO_DZ dv_do_dz;
  typename MODELS::OMEGA_Z omega_z;

  y3_cluster::IntegrationRange lnM_ir_;
  std::array<y3_cluster::IntegrationRange, NRICHNESS> lo_ir_;
  y3_cluster::IntegrationRange lt_ir_;
  y3_cluster::IntegrationRange lc_ir_;
  std::array<y3_cluster::IntegrationRange, NREDSHIFT> zo_ir_;
  y3_cluster::IntegrationRange zt_ir_;
  y3_cluster::IntegrationRange R_ir_;
  y3_cluster::IntegrationRange A_ir_;
  y3_cluster::IntegrationRange theta_ir_;

  std::array <double, NRADII> r;
public:
  // A Gamma_T_Integrand object is constructed by passing in the bunch of
  // callable objects (function pointers or callable class instances)  that
  // specify the various terms of the integrand.
  Gamma_T_Integrand(double fcen,
                    double msci,
                    typename MODELS::MOR mor,
                    typename MODELS::LO_LC lo_lc,
                    typename MODELS::LC_LT lc_lt,
                    typename MODELS::ZO_ZT zo_zt,
                    typename MODELS::ROFFSET roffset,
                    typename MODELS::T_CEN T_cen,
                    typename MODELS::T_MIS T_mis,
                    typename MODELS::A_CEN A_cen,
                    typename MODELS::A_MIS A_mis,
                    typename MODELS::HMF hmf,
                    typename MODELS::DEL_SIG del_sig,
                    typename MODELS::DV_DO_DZ dv_do_dz,
                    typename MODELS::OMEGA_Z omega_z,
                    y3_cluster::IntegrationRange lnM_ir,
                    std::array<y3_cluster::IntegrationRange, NRICHNESS> lo_ir,
                    y3_cluster::IntegrationRange lt_ir,
                    y3_cluster::IntegrationRange lc_ir,
                    std::array<y3_cluster::IntegrationRange, NREDSHIFT> zo_ir,
                    y3_cluster::IntegrationRange zt_ir,
                    y3_cluster::IntegrationRange R_ir,
                    y3_cluster::IntegrationRange A_ir,
                    y3_cluster::IntegrationRange theta_ir,
                    std::array<double, NRADII> const& rarray)
    : fcen_(fcen)
    , msci_(msci)
    , mor(mor)
    , lo_lc(lo_lc)
    , lc_lt(lc_lt)
    , zo_zt(zo_zt)
    , roffset(roffset)
    , T_cen(T_cen)
    , T_mis(T_mis)
    , A_cen(A_cen)
    , A_mis(A_mis)
    , hmf(hmf)
    , del_sig(del_sig)
    , dv_do_dz(dv_do_dz)
    , omega_z(omega_z)
    , lnM_ir_(lnM_ir)
    , lo_ir_(lo_ir)
    , lt_ir_(lt_ir)
    , lc_ir_(lc_ir)
    , zo_ir_(zo_ir)
    , zt_ir_(zt_ir)
    , R_ir_(R_ir)
    , A_ir_(A_ir)
    , theta_ir_(theta_ir)	
    , r(rarray)
  {}

  template <std::size_t NRanges>
  static inline std::array<y3_cluster::IntegrationRange, NRanges>
  _get_ranges(cosmosis::DataBlock db, std::string name)
  {
      std::array<std::size_t, NRanges> r;
      for (auto i = 0u; i < NRanges; i++)
          r[i] = i;
      auto mins = db.view<std::vector<double>>("integration_ranges", name + "_min");
      auto maxs = db.view<std::vector<double>>("integration_ranges", name + "_min");
      return y3_cluster::transform(r, [&](std::size_t i) {
                return y3_cluster::IntegrationRange{mins[i], maxs[i]};
              });
  }

  // Alternatively, can automatically construct each model component given a datablock.
  Gamma_T_Integrand(cosmosis::DataBlock& sample)
    : mor(sample)
    , lo_lc(sample)
    , lc_lt(sample)
    , zo_zt(sample)
    , roffset(sample)
    , T_cen(sample)
    , T_mis(sample)
    , A_cen(sample)
    , A_mis(sample)
    , hmf(sample)
    , del_sig(sample)
    , dv_do_dz(sample)
    , omega_z(sample)
    , lnM_ir_(sample, "lnM")
    , lo_ir_(_get_ranges<NRICHNESS>(sample, "lo"))
    , lt_ir_(sample, "lt")
    , lc_ir_(sample, "lc")
    , zo_ir_(_get_ranges<NREDSHIFT>(sample, "zo"))
    , zt_ir_(sample, "zt")
    , R_ir_(sample, "R")
    , A_ir_(sample, "A")
    , theta_ir_(sample, "theta")
    //, r(sample)
  // TODO: Switch this out for a general method soon!
  {
    for (std::size_t i = 0; i < 10; i++) {
      r[i] = 5.0 * (i + 0.01);
    }
  }

  template<std::size_t NEW_NRICHNESS, std::size_t NEW_NREDSHIFT>
  Gamma_T_Integrand<MODELS, NRADII, NEW_NRICHNESS, NEW_NREDSHIFT>
  with_bins(std::array<y3_cluster::IntegrationRange, NEW_NRICHNESS> new_lir,
            std::array<y3_cluster::IntegrationRange, NEW_NREDSHIFT> new_zir)
  {
      return {fcen_,
              msci_,
              mor,
              lo_lc,
              lc_lt,
              zo_zt,
              roffset,
              T_cen,
              T_mis,
              A_cen,
              A_mis,
              hmf,
              del_sig,
              dv_do_dz,
              omega_z,
              lnM_ir_,
              // Different lir
              new_lir,
              lt_ir_,
              lc_ir_,
              // Different zir
              new_zir,
              zt_ir_,
              R_ir_,
              A_ir_,
              theta_ir_,
              r};
  }

  template<typename Fjn, typename Fjg, typename Fnm, typename Fgr>
  std::array<double, NRICHNESS * NREDSHIFT * (NRADII+2)>
  integrand_common(double lt,
                   double zt,
                   double lnM,
                   // Jacobian for N term
                   Fjn jacob_N,
                   // Jacobian for Gamma term
                   Fjg jacob_G,
                   Fnm N_multiplier,
                   // Radially dependent function
                   Fgr gamma_radial_dep) const
  {
    std::array<double, NRICHNESS * NREDSHIFT * (NRADII+2)> return_arr;

    auto const hmf_v = hmf(lnM, zt);
    auto const mor_v = mor(lt, lnM, zt);
    auto const dv_do_dz_v = dv_do_dz(zt);
    auto const omega_z_v = omega_z(zt);

    for (std::size_t loi = 0; loi < NRICHNESS; loi++) {
        auto const richness_bin_start = loi * NREDSHIFT * (NRADII + 2);

        for (std::size_t zoi = 0; zoi < NREDSHIFT; zoi++) {
            // Zo does not actually need to be integrated over
            double const zomin = zo_ir_[zoi].transform(0.0);
            double const zomax = zo_ir_[zoi].transform(1.0);
            auto const zo_zt_v = zo_zt(zomin, zomax, zt);

            // These will eventually be passed by CosmoSIS
            double m_shear = 0.0;
            double sig_crit = 1.0;
            // This is the lambda-redshift bin weight that we don't fully understand
            double w = 1.0;

            // eq. (25)
            double const N_int = omega_z_v * dv_do_dz_v * zo_zt_v * hmf_v * mor_v;
            double const N_mult = N_multiplier(loi);
            // eq. (24)
            double const N = jacob_N(loi) * N_int * N_mult;
            double const Nw = N * w;

            // eq. (29)
            auto const gamma_t_int = jacob_G(loi) * N_int * w;

            // eq. (28)
            auto const gamma_t = y3_cluster::transform(r,
	                   [m_shear, sig_crit, gamma_t_int, N_mult, gamma_radial_dep]
                           (double radius) {
                               // Nw intentionally left out - returned in return_arr to be used further on
                               return (1.0 + m_shear) / sig_crit
                                       * gamma_t_int * N_mult * gamma_radial_dep(radius);
                           });

            auto redshift_bin_start = richness_bin_start + zoi * (NRADII + 2);

            std::copy_n( gamma_t.begin(), gamma_t.size(), &return_arr[redshift_bin_start] );
            return_arr[redshift_bin_start + NRADII] = N;
            return_arr[redshift_bin_start + NRADII + 1] = Nw;
        }
    }

    return return_arr;
  }

  /* Miscentered part of integration
   * Term dictionary -
   * * lo - \lambda^{obs}
   * * lc - \lambda^{cen}
   * * lt - \lambda^{true}
   * * zo - z^{obs}
   * * zt - z^{true}
   * * R - R
   * * lnM - ln(M)
   * * A - ???
   * */
  std::array<double, NRICHNESS * NREDSHIFT * (NRADII+2)>
  miscentered(double scaled_lo,
              double scaled_lc,
              double scaled_lt,
              double scaled_zt,
              double scaled_R,
              double scaled_lnM,
              double scaled_A,
              double scaled_theta) const
  {
    // We probably should factor out the common subexpressions, rather than
    // relying upon the optimizer to do a perfect job of this for us. This
    // seems to be the intent of the commented-out code below.
    auto const lnM = lnM_ir_.transform(scaled_lnM);
    auto const lt = lt_ir_.transform(scaled_lt);
    auto const lc = lc_ir_.transform(scaled_lc);
    auto const zt = zt_ir_.transform(scaled_zt);
    auto const R = R_ir_.transform(scaled_R);
    auto const A = A_ir_.transform(scaled_A);
    auto const theta = theta_ir_.transform(scaled_theta);

    auto jacob_N = [=](std::size_t loi) {
        return lnM_ir_.jacobian() * lo_ir_[loi].jacobian()
               * lt_ir_.jacobian() * lc_ir_.jacobian()
               * zt_ir_.jacobian()
               * R_ir_.jacobian();
    };
    auto jacob_G = [=](std::size_t loi) {
        return lnM_ir_.jacobian() * lo_ir_[loi].jacobian()
               * lt_ir_.jacobian() * lc_ir_.jacobian()
               * zt_ir_.jacobian()
               * R_ir_.jacobian() * A_ir_.jacobian()
               * theta_ir_.jacobian();
    };

    // eq. (27)
    auto N_mis = [=](std::size_t loi) {
        auto const lo = lo_ir_[loi].transform(scaled_lo);
        return (1.0 - fcen_) * lo_lc(lo, lc, R) * lc_lt(lc, lt, zt) * roffset(R);
    };

    // eq. (30)
    // For the following lambda function, `radius` corresponds to what is called
    // `R` in the paper, and `R` corresponds to what is called `R_{mis}` in the
    // paper
    // eq. (31)
    auto gamma_t_mis = [this, N_mis, A, lnM, R, theta, zt](double radius) {
        double const adjusted_R = std::sqrt(radius*radius + R*R + 2*R*radius * std::cos(theta));
        return (1.0 / 6.28318530718)
               * exp(A * T_cen(adjusted_R, lnM))/A_ir_.jacobian()
               * del_sig(adjusted_R, lnM, zt);
    };

    return integrand_common(lt,
                            zt,
                            lnM,
                            jacob_N,
                            jacob_G,
                            N_mis,
                            gamma_t_mis);
  }

  /* Centered part of integration
   * Term dictionary -
   * * lo - \lambda^{obs}
   * * lc - \lambda^{cen}
   * * lt - \lambda^{true}
   * * zo - z^{obs}
   * * zt - z^{true}
   * * R - R
   * * lnM - ln(M)
   * * A - ???
   * */
  std::array<double, NRICHNESS * NREDSHIFT * (NRADII+2)>
  centered(double scaled_lo,
           double scaled_lt,
           double scaled_zt,
           double scaled_lnM,
           double scaled_A) const
  {
    // Necessary terms
    auto const lnM = lnM_ir_.transform(scaled_lnM);
    auto const lt = lt_ir_.transform(scaled_lt);
    auto const zt = zt_ir_.transform(scaled_zt);
    auto const A = A_ir_.transform(scaled_A);

    auto jacob_N = [=](std::size_t loi) {
        return lnM_ir_.jacobian() * lo_ir_[loi].jacobian()
               * lt_ir_.jacobian()
               * zt_ir_.jacobian();
    };

    auto jacob_G = [=](std::size_t loi) {
        return lnM_ir_.jacobian() * lo_ir_[loi].jacobian()
               * lt_ir_.jacobian()
               * zt_ir_.jacobian()
               * A_ir_.jacobian();
    };

    // eq. (26)
    auto N_cen = [=](std::size_t loi) {
        auto const lo = lo_ir_[loi].transform(scaled_lo);
        return lc_lt(lo, lt, zt) * fcen_;
    };

    // eq. (30)
    // For the following lambda function, `radius` corresponds to what is called
    // `R` in the paper
    auto gamma_t_cen = [this, N_cen, A, lnM, zt](double radius) {
        return exp(A * T_cen(radius, lnM)) / A_ir_.jacobian() * del_sig(radius, lnM, zt);
    };

    return integrand_common(lt,
                            zt,
                            lnM,
                            jacob_N,
                            jacob_G,
                            N_cen,
                            gamma_t_cen);
  }

  template<typename Integrator>
  cubacpp::integration_results<NRADII + 2>
  integrate_centered(Integrator i, double epsrel, double epsabs)
  {
      return i.integrate([this](double scaled_lo, double scaled_lt, double scaled_zt,
                                double scaled_lnM, double scaled_A) {
                            return centered(scaled_lo, scaled_lt, scaled_zt, scaled_lnM, scaled_A);
                         },
                         epsrel, epsabs);
  }

  template<typename Integrator>
  cubacpp::integration_results<NRADII + 2>
  integrate_miscentered(Integrator i, double epsrel, double epsabs)
  {
      return i.integrate([this](double scaled_lo, double scaled_lc, double scaled_lt,
                                double scaled_zt, double scaled_R, double scaled_lnM,
                                double scaled_A, double scaled_theta) {
                            return miscentered(scaled_lo, scaled_lc, scaled_lt,
                                               scaled_zt, scaled_R, scaled_lnM,
                                               scaled_A, scaled_theta);
                         },
                         epsrel, epsabs);
  }
};

template <typename MODELS, std::size_t NRADII=10, std::size_t NRICHNESS=1, std::size_t NREDSHIFT=1>
Gamma_T_Integrand<MODELS, NRADII, NRICHNESS, NREDSHIFT>
make_gamma_t_integrand(double fcen,
                       double msci,
                       typename MODELS::MOR mor,
                       typename MODELS::LO_LC lo_lc,
                       typename MODELS::LC_LT lc_lt,
                       typename MODELS::ZO_ZT zo_zt,
                       typename MODELS::ROFFSET roffset,
                       typename MODELS::T_CEN t_cen,
                       typename MODELS::T_MIS t_mis,
                       typename MODELS::A_CEN a_cen,
                       typename MODELS::A_MIS a_mis,
                       typename MODELS::HMF hmf,
                       typename MODELS::DEL_SIG del_sig,
                       typename MODELS::DV_DO_DZ dv_do_dz,
                       typename MODELS::OMEGA_Z omega_z,
                       std::array<y3_cluster::IntegrationRange, NRICHNESS> lo_ir,
                       std::array<y3_cluster::IntegrationRange, NREDSHIFT> zo_ir)
{
   y3_cluster::IntegrationRange lnM_ir{std::log(6.e12), std::log(4.e16)};
   y3_cluster::IntegrationRange lt_ir{2.0, 100}; // we should adjust lt, lc and lnM ranges according to the bin
   y3_cluster::IntegrationRange lc_ir{2.0, 100};
   y3_cluster::IntegrationRange zt_ir{0.05, 0.3};
   y3_cluster::IntegrationRange R_ir{0., 2.0};
   y3_cluster::IntegrationRange A_ir{-0.01, 0.01};
   y3_cluster::IntegrationRange theta_ir{0.,6.28318530718};

   std::array<double, NRADII> rarray; 
   for (std::size_t i = 0; i < NRADII; i++)
       rarray[i] = 0.1 * (i + 0.1);

   return {fcen,
           msci,
           mor,
           lo_lc,
           lc_lt,
           zo_zt,
           roffset,
           t_cen,
           t_mis,
           a_cen,
           a_mis,
           hmf,
           del_sig,
           dv_do_dz,
           omega_z,
           lnM_ir,
           lo_ir,
           lt_ir,
           lc_ir,
           zo_ir,
           zt_ir,
           R_ir,
           A_ir,
           theta_ir,
           rarray };
}

template <std::size_t NRADII>
struct Gamma_T_Integrated_Bin_Result {
    y3_cluster::IntegrationRange lo_ir, zo_ir;
    std::array<double, NRADII> radius;
    std::array<double, NRADII> gamma_ts;
    std::array<double, NRADII> gamma_t_errors;
    std::array<double, NRADII> gamma_t_probs;
    double N, N_error, N_prob,
           Nw, Nw_error, Nw_prob;

    Gamma_T_Integrated_Bin_Result() : lo_ir{0.0, 1.0}, zo_ir{0.0, 1.0} {}

    template<typename MODELS, std::size_t NRICHNESS, std::size_t NREDSHIFT>
    Gamma_T_Integrated_Bin_Result(std::size_t which_richness,
                                  std::size_t which_redshift,
                                  const Gamma_T_Integrand<MODELS, NRADII, NRICHNESS, NREDSHIFT>& gt,
                                  const cubacpp::integration_results<NRICHNESS * NREDSHIFT * (NRADII + 2)>& results)
        : lo_ir(gt.lo_ir_[which_richness])
        , zo_ir(gt.zo_ir_[which_redshift])
        , radius(gt.r)
    {
        const auto base = (which_richness * NREDSHIFT + which_redshift) * (NRADII + 2);

        for (auto i = 0u; i < NRADII; i++) {
            gamma_ts[i] = results.value[base + i];
            gamma_t_errors[i] = results.error[base + i];
            gamma_t_probs[i] = results.prob[base + i];
        }

        N = results.value[base + NRADII];
        Nw = results.value[base + NRADII + 1];

        N_error = results.error[base + NRADII];
        Nw_error = results.error[base + NRADII + 1];

        N_prob = results.prob[base + NRADII];
        Nw_prob = results.prob[base + NRADII + 1];
    }
};

template<typename MODELS, std::size_t NRADII, std::size_t NRICHNESS = 1, std::size_t NREDSHIFT = 1>
std::array<Gamma_T_Integrated_Bin_Result<NRADII>, NRICHNESS * NREDSHIFT>
make_gamma_t_integrated_bins(const Gamma_T_Integrand<MODELS, NRADII, NRICHNESS, NREDSHIFT>& gt,
                             const cubacpp::integration_results<NRICHNESS * NREDSHIFT * (NRADII + 2)>& results)
{
    std::array<Gamma_T_Integrated_Bin_Result<NRADII>, NRICHNESS * NREDSHIFT> return_arr;

    for (auto loi = 0u; loi < NRICHNESS; loi++) {
        for (auto zoi = 0u; zoi < NREDSHIFT; zoi++) {
            return_arr[loi * NREDSHIFT + zoi] = Gamma_T_Integrated_Bin_Result<NRADII>(loi, zoi, gt, results);
        }
    }

    return return_arr;
}

} // namespace y3_cluster

#endif
