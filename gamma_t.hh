#ifndef Y3_CLUSTER_GAMMA_T_HH
#define Y3_CLUSTER_GAMMA_T_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "integration_range.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <test/transform.hh>

// This class template is based
// on https://www.overleaf.com/13697016cyvvqqfchfbg#/52989522/, and the example
// provided by Spencer Everett.

template <typename MODELS>
class Gamma_T_Integrand {
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
  typename MODELS::DEL_SIG_MIS del_sig_mis;
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

  // How to declare this if reading NRADII from datablock?
  static const std::size_t NRADII = 10;
  std::array<double, NRADII> r;

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
                    typename MODELS::DEL_SIG_MIS del_sig_mis,
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
    , del_sig_mis(del_sig_mis)
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
    , r(rarray)
  {}

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
    , del_sig_cen(sample)
    , del_sig_mis(sample)
    , dv_do_dz(sample)
    , omega_z(sample)
    , lnM_ir_(sample)
    , lo_ir_(sample)
    , lt_ir_(sample)
    , lc_ir_(sample)
    , zo_ir_(sample)
    , zt_ir_(sample)
    , R_ir_(sample)
    , A_ir_(sample)
    , r(sample)
  {}

  // The function call operator -- this is the function to be integrated.
  std::array<double, NRADII + 2>
  operator()(double scaled_lo,
             double scaled_lc,
             double scaled_lt,
             double scaled_zo,
             double scaled_zt,
             double scaled_R,
             double scaled_lnM,
             double scaled_A) const
  {
    // We probably should factor out the common subexpressions, rather than
    // relying upon the optimizer to do a perfect job of this for us. This
    // seems to be the intent of the commented-out code below.
    using std::exp;
    auto const lnM = lnM_ir_.transform(scaled_lnM);
    auto const lo = lo_ir_.transform(scaled_lo);
    auto const lt = lt_ir_.transform(scaled_lt);
    auto const lc = lc_ir_.transform(scaled_lc);
    auto const zo = zo_ir_.transform(scaled_zo);
    auto const zt = zt_ir_.transform(scaled_zt);
    auto const R = R_ir_.transform(scaled_R);
    auto const A = A_ir_.transform(scaled_A);

    auto const hmf_v = hmf(lnM, zt);
    auto const zo_zt_v = zo_zt(zo, zt);
    auto const lc_lt_v = lc_lt(lc, lt, zt);
    auto const mor_v = mor(lt, lnM, zt);
    auto const dv_do_dz_v = dv_do_dz(zt);
    auto const omega_z_v = omega_z(zt);

    // These will eventually be passed by CosmoSIS
    double m_shear = 1.0;
    double sig_crit_inv = 1.0;
    // This is the lambda-redshift bin weight that we don't fully understand
    double w = 1.0;

    // The evaluation below follows the convention set in main overleaf document
    // The evaluation is for Y3 likelihood
    // puttign together the return vector
    double const jacob = lnM_ir_.jacobian() * lo_ir_.jacobian() *
                         lt_ir_.jacobian() * lc_ir_.jacobian() *
                         zo_ir_.jacobian() * zt_ir_.jacobian() *
                         R_ir_.jacobian() * A_ir_.jacobian();
    double const N = jacob * omega_z_v * dv_do_dz_v * zo_zt_v * hmf_v * mor_v *
                     lc_lt_v *
                     (fcen_ + (1.0 - fcen_) * roffset(R) * lo_lc(lo, lc, R));
    double const Nw = jacob * N * w;

    auto const gamma_t_int =
      jacob * omega_z_v * dv_do_dz_v * zo_zt_v * hmf_v * mor_v * w * lc_lt_v;
    auto gamma_t_cen = [this, lnM, zt, A](double radius) {
      return fcen_ * del_sig_cen(radius, lnM, zt) * exp(A * T_cen(radius, lnM));
    };
    auto gamma_t_mis = [this, lnM, A, R, lo, lc](double radius) {
      return (1.0 - fcen_) * roffset(R) * lo_lc(lo, lc, R) *
             del_sig_mis(radius, lnM, R) * exp(A * T_cen(radius, lnM));
    };
    auto const gamma_t = y3_cluster::transform(
      r,
      [gamma_t_cen, gamma_t_mis, m_shear, sig_crit_inv, gamma_t_int](
        double radius) {
        return (1.0 + m_shear) / sig_crit_inv * gamma_t_int *
               (gamma_t_cen(radius) + gamma_t_mis(radius));
      });

    std::array<double, NRADII + 2> return_arr;
    std::copy_n(gamma_t.begin(), gamma_t.size(), return_arr.begin());
    return_arr[NRADII] = N;
    return_arr[NRADII + 1] = Nw;

    return return_arr;
  }
};

// Helper function to create a 'standard' gamma_T integrand
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
                       typename MODELS::DEL_SIG_MIS del_sig_mis,
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
  std::array<double, 10> rarray; // can I pass a vector here?
  for (std::size_t i = 0; i < 10; i++) {
    rarray[i] = 5.0 * (i + 0.01);
  }
  return {fcen,     msci,    mor,    lo_lc, lc_lt, zo_zt,       roffset,
          t_cen,    t_mis,   a_cen,  a_mis, hmf,   del_sig_cen, del_sig_mis,
          dv_do_dz, omega_z, lnM_ir, lo_ir, lt_ir, lc_ir,       zo_ir,
          zt_ir,    R_ir,    A_ir,   rarray};
}

template <typename MODELS>
Gamma_T_Integrand<MODELS>
make_gamma_t_integrand(cosmosis::DataBlock& sample)
{
  // TODO: Allow for integration ranges from the datablock!
  // Defaults?
  y3_cluster::IntegrationRange lnM_ir{std::log(5.e11), std::log(1.e17)};
  y3_cluster::IntegrationRange lt_ir{1.0, 100};
  y3_cluster::IntegrationRange lc_ir{1.0, 100};
  y3_cluster::IntegrationRange zt_ir{0.1, 0.3};
  y3_cluster::IntegrationRange R_ir{0., 1.0};
  y3_cluster::IntegrationRange A_ir{-1.0, 1.0};
  std::array<double, 10> rarray; // can I pass a vector here?
  for (std::size_t i = 0; i < 10; i++) {
    rarray[i] = 5.0 * (i + 0.01);
  }
  return {fcen,     msci,    mor,    lo_lc, lc_lt, zo_zt,       roffset,
          t_cen,    t_mis,   a_cen,  a_mis, hmf,   del_sig_cen, del_sig_mis,
          dv_do_dz, omega_z, lnM_ir, lo_ir, lt_ir, lc_ir,       zo_ir,
          zt_ir,    R_ir,    A_ir,   rarray};
}

#endif
