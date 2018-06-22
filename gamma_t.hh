#ifndef Y3_CLUSTER_GAMMA_T_HH
#define Y3_CLUSTER_GAMMA_T_HH

#include "integration_range.hh"

#include <array>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <test/transform.hh>
// This class template is based
// on https://www.overleaf.com/13697016cyvvqqfchfbg#/52989522/, and the example
// provided by Spencer Everett.



template<typename MOR_,
         typename LO_LC_,
         typename LC_LT_,
         typename ZO_ZT_,
         typename ROFFSET_,
         typename T_CEN_,
         typename T_MIS_,
         typename A_CEN_,
         typename A_MIS_,
         typename HMF_,
         typename DEL_SIG_CEN_,
         typename DV_DO_DZ_,
         typename OMEGA_Z_>
struct Models {
    using MOR = MOR_;
    using LO_LC = LO_LC_;
    using LC_LT = LC_LT_;
    using ZO_ZT = ZO_ZT_;
    using ROFFSET = ROFFSET_;
    using T_CEN = T_CEN_;
    using T_MIS = T_MIS_;
    using A_CEN = A_CEN_;
    using A_MIS = A_MIS_;
    using HMF = HMF_;
    using DEL_SIG_CEN = DEL_SIG_CEN_;
    using DV_DO_DZ = DV_DO_DZ_;
    using OMEGA_Z = OMEGA_Z_;
};

template <typename MODELS>
class Gamma_T_Integrand {
public:
#ifndef OVERRIDE_NRADII
  static const std::size_t NRADII = 10;
#else
  static const std::size_t NRADII = OVERRIDE_NRADII;
#endif
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
  typename MODELS::DEL_SIG_CEN del_sig_cen;
  typename MODELS::DV_DO_DZ dv_do_dz;
  typename MODELS::OMEGA_Z omega_z;

  y3_cluster::IntegrationRange lnM_ir_;
  y3_cluster::IntegrationRange lo_ir_;
  y3_cluster::IntegrationRange lt_ir_;
  y3_cluster::IntegrationRange lc_ir_;
  y3_cluster::IntegrationRange zo_ir_;
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
                    typename MODELS::DEL_SIG_CEN del_sig_cen,
                    typename MODELS::DV_DO_DZ dv_do_dz,
                    typename MODELS::OMEGA_Z omega_z,
                    y3_cluster::IntegrationRange lnM_ir,
                    y3_cluster::IntegrationRange lo_ir,
                    y3_cluster::IntegrationRange lt_ir,
                    y3_cluster::IntegrationRange lc_ir,
                    y3_cluster::IntegrationRange zo_ir,
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
    , del_sig_cen(del_sig_cen)
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

  template<typename F>
  std::array<double, NRADII+2>
  integrand_common(double lt,
                   // double zo,
                   double zt,
                   double lnM,
                   // Jacobian for N term
                   double jacob_N,
                   // Jacobian for Gamma term
                   double jacob_G,
                   double N_multiplier,
                   // Radially dependent function
                   F gamma_radial_dep
                   ) const
  {
    // Zo does not actually need to be integrated over
    double const zomin = zo_ir_.transform(0.0);
    double const zomax = zo_ir_.transform(1.0);

    auto const hmf_v = hmf(lnM, zt);
    auto const zo_zt_v = zo_zt(zomin, zomax, zt);
    auto const mor_v = mor(lt, lnM, zt);
    auto const dv_do_dz_v = dv_do_dz(zt);
    auto const omega_z_v = omega_z(zt);


    // These will eventually be passed by CosmoSIS
    double m_shear = 1.0;
    double sig_crit = 1.0;
    // This is the lambda-redshift bin weight that we don't fully understand
    double w = 1.0;

    // eq. (25)
    double const N_int = omega_z_v * dv_do_dz_v * zo_zt_v * hmf_v * mor_v;
    // eq. (24)
    double const N = jacob_N * N_int * N_multiplier;
    double const Nw = N * w;

    // eq. (29)
    auto const gamma_t_int = jacob_G * N_int * w;

    // eq. (28)
    auto const gamma_t = y3_cluster::transform(r,
	           [m_shear, sig_crit, gamma_t_int, gamma_radial_dep]
                   (double radius) {
                       // Nw intentionally left out - returned in return_arr to be used further on
                       return (1.0 + m_shear) / sig_crit
                               * gamma_t_int * gamma_radial_dep(radius);
                   });

    std::array<double, NRADII+2> return_arr;
    std::copy_n( gamma_t.begin(), gamma_t.size(),  return_arr.begin() );
    return_arr[NRADII] = N;
    return_arr[NRADII+1] = Nw;

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
  std::array<double, NRADII+2>
  miscentered(double scaled_lo,
              double scaled_lc,
              double scaled_lt,
              // double scaled_zo,
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
    auto const lo = lo_ir_.transform(scaled_lo);
    auto const lt = lt_ir_.transform(scaled_lt);
    auto const lc = lc_ir_.transform(scaled_lc);
    auto const zt = zt_ir_.transform(scaled_zt);
    auto const R = R_ir_.transform(scaled_R);
    auto const A = A_ir_.transform(scaled_A);
    auto const theta = theta_ir_.transform(scaled_theta);

    double const jacob_N = lnM_ir_.jacobian() * lo_ir_.jacobian()
                         * lt_ir_.jacobian() * lc_ir_.jacobian()
                         //* zo_ir_.jacobian()
                         * zt_ir_.jacobian()
                         * R_ir_.jacobian();
    double const jacob_G = lnM_ir_.jacobian() * lo_ir_.jacobian()
                         * lt_ir_.jacobian() * lc_ir_.jacobian()
                         //* zo_ir_.jacobian()
                         * zt_ir_.jacobian()
                         * R_ir_.jacobian() * A_ir_.jacobian()
                         * theta_ir_.jacobian();

    // eq. (27)
    double const N_mis = (1.0 - fcen_) * lo_lc(lo, lc, R) * lc_lt(lc, lt, zt) * roffset(R);

    // eq. (30)
    // For the following lambda function, `radius` corresponds to what is called
    // `R` in the paper, and `R` corresponds to what is called `R_{mis}` in the
    // paper
    // eq. (31)
    auto gamma_t_mis = [this, N_mis, A, lnM, R, theta, zt](double radius) {
        double const adjusted_R = std::sqrt(radius*radius + R*R + 2*R*radius * std::cos(theta));
        return (N_mis / 6.28318530718)
               * exp(A * T_cen(adjusted_R, lnM))
               * del_sig_cen(adjusted_R, lnM, zt);
    };

    return integrand_common(lt,
                            // zo_ir_.transform(scaled_zo),
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
  std::array<double, NRADII+2>
  centered(double scaled_lo,
           double scaled_lt,
           // double scaled_zo,
           double scaled_zt,
           double scaled_lnM,
           double scaled_A) const
  {
    // Necessary terms
    auto const lnM = lnM_ir_.transform(scaled_lnM);
    auto const lo = lo_ir_.transform(scaled_lo);
    auto const lt = lt_ir_.transform(scaled_lt);
    auto const zt = zt_ir_.transform(scaled_zt);
    auto const A = A_ir_.transform(scaled_A);

    double const jacob_N = lnM_ir_.jacobian() * lo_ir_.jacobian()
                         * lt_ir_.jacobian()
                         //* zo_ir_.jacobian()
                         * zt_ir_.jacobian();
    double const jacob_G = lnM_ir_.jacobian() * lo_ir_.jacobian()
                         * lt_ir_.jacobian()
                         //* zo_ir_.jacobian()
                         * zt_ir_.jacobian()
                         * A_ir_.jacobian();

    // eq. (26)
    auto const lo_lt_v = lc_lt(lo, lt, zt);
    double const N_cen = lo_lt_v * fcen_;

    // eq. (30)
    // For the following lambda function, `radius` corresponds to what is called
    // `R` in the paper
    auto gamma_t_cen = [this, N_cen, A, lnM, zt](double radius) {
        return N_cen * exp(A * T_cen(radius, lnM)) * del_sig_cen(radius, lnM, zt);
    };

    return integrand_common(lt,
                            // zo_ir_.transform(scaled_zo),
                            zt,
                            lnM,
                            jacob_N,
                            jacob_G,
                            N_cen,
                            gamma_t_cen);
  }
};

template <typename MODELS>
Gamma_T_Integrand<MODELS>
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
                       typename MODELS::DEL_SIG_CEN del_sig_cen,
                       typename MODELS::DV_DO_DZ dv_do_dz,
                       typename MODELS::OMEGA_Z omega_z,
                       y3_cluster::IntegrationRange lo_ir,
                       y3_cluster::IntegrationRange zo_ir)
{
   y3_cluster::IntegrationRange lnM_ir{std::log(5.e11), std::log(1.e17)};
   y3_cluster::IntegrationRange lt_ir{1.0, 100};
   y3_cluster::IntegrationRange lc_ir{1.0, 100};
   y3_cluster::IntegrationRange zt_ir{0.1, 0.3};
   y3_cluster::IntegrationRange R_ir{0., 1.0};
   y3_cluster::IntegrationRange A_ir{-1.0, 1.0};
   y3_cluster::IntegrationRange theta_ir{0.,6.28318530718};

   std::size_t const NRADII = Gamma_T_Integrand<MODELS>::NRADII;
   std::array<double, NRADII> rarray; // can I pass a vector here?
   for ( std::size_t i = 0; i < NRADII; i++ ) {rarray[ i ] = 0.1*(i+0.1);}
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
           del_sig_cen,
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

#endif
