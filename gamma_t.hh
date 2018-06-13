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

template <typename MOR,
          typename LO_LC,
          typename LC_LT,
          typename ZO_ZT,
          typename ROFFSET,
          typename T_CEN,
          typename T_MIS,
          typename A_CEN,
          typename A_MIS,
          typename HMF,
          typename DEL_SIG_CEN,
          /* typename DEL_SIG_MIS,*/
          typename DV_DO_DZ,
          typename OMEGA_Z>
class Gamma_T_Integrand {
private:
  double fcen_;
  double msci_;

  MOR mor;
  LO_LC lo_lc;
  LC_LT lc_lt;
  ZO_ZT zo_zt;
  ROFFSET roffset;
  T_CEN T_cen;
  T_MIS T_mis;
  A_CEN A_cen;
  A_MIS A_mis;
  HMF hmf;
  DEL_SIG_CEN del_sig_cen;
  /*DEL_SIG_MIS del_sig_mis;*/
  DV_DO_DZ dv_do_dz;
  OMEGA_Z omega_z;

  y3_cluster::IntegrationRange lnM_ir_;
  y3_cluster::IntegrationRange lo_ir_;
  y3_cluster::IntegrationRange lt_ir_;
  y3_cluster::IntegrationRange lc_ir_;
  y3_cluster::IntegrationRange zo_ir_;
  y3_cluster::IntegrationRange zt_ir_;
  y3_cluster::IntegrationRange R_ir_;
  y3_cluster::IntegrationRange A_ir_;
  y3_cluster::IntegrationRange theta_ir_; /* */

  static const std::size_t NRADII = 10;
  std::array <double, NRADII> r;
public:
  // A Gamma_T_Integrand object is constructed by passing in the bunch of
  // callable objects (function pointers or callable class instances)  that
  // specify the various terms of the integrand.
  Gamma_T_Integrand(double fcen,
                    double msci,
                    MOR mor,
                    LO_LC lo_lc,
                    LC_LT lc_lt,
                    ZO_ZT zo_zt,
                    ROFFSET roffset,
                    T_CEN T_cen,
                    T_MIS T_mis,
                    A_CEN A_cen,
                    A_MIS A_mis,
                    HMF hmf,
                    DEL_SIG_CEN del_sig_cen,
                    /*DEL_SIG_MIS del_sig_mis,*/
                    DV_DO_DZ dv_do_dz,
                    OMEGA_Z omega_z,
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
    /*, del_sig_mis(del_sig_mis)*/
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

  // The function call operator -- this is the function to be integrated.
  /* Term dictionary -
   * * lo - \lambda^{obs}
   * * lc - \lambda^{cen}
   * * lt - \lambda^{true}
   * * zo - z^{obs}
   * * zt - z^{true}
   * * R - R
   * * lnM - ln(M)   // (why log???)
   * * A - ???
   * */
  std::array<double, NRADII+2>
  operator()(double scaled_lo,
             double scaled_lc,
             double scaled_lt,
             double scaled_zo,
             double scaled_zt,
             double scaled_R,
             double scaled_lnM,
             double scaled_A,
	     double scaled_theta) const
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
    auto const theta = theta_ir_.transform(scaled_theta);

    auto const hmf_v = hmf(lnM, zt);
    auto const zo_zt_v = zo_zt(zo, zt);
    auto const lc_lt_v = lc_lt(lc, lt, zt);
    auto const lo_lt_v = lc_lt(lo, lt, zt);
    auto const mor_v = mor(lt, lnM, zt);
    auto const dv_do_dz_v = dv_do_dz(zt);
    auto const omega_z_v = omega_z(zt);


    // These will eventually be passed by CosmoSIS
    double m_shear = 1.0;
    double sig_crit = 1.0;
    // This is the lambda-redshift bin weight that we don't fully understand
    double w = 1.0;

    // The evaluation below follows the convention set in main overleaf document
    //The evaluation is for Y3 likelihood
    // putting together the return vector
    double const lc_jacob = lc_ir_.jacobian();
    double const jacob_N = lnM_ir_.jacobian() * lo_ir_.jacobian()
                        * lt_ir_.jacobian() * lc_ir_.jacobian()
                        * zo_ir_.jacobian() * zt_ir_.jacobian()
                        * R_ir_.jacobian();
    double const jacob = lnM_ir_.jacobian() * lo_ir_.jacobian()
                       * lt_ir_.jacobian() * lc_ir_.jacobian()
                       * zo_ir_.jacobian() * zt_ir_.jacobian()
                       * R_ir_.jacobian() * A_ir_.jacobian()
                       * theta_ir_.jacobian();

    /* eq. (25) */
    double const N = jacob_N * omega_z_v * dv_do_dz_v * zo_zt_v * hmf_v * mor_v
                   /* eq. (26) + eq. (27) */
                   * (lo_lt_v * fcen_ / lc_jacob + lc_lt_v * (1.0 - fcen_) * roffset(R) * lo_lc(lo, lc, R));
    double const Nw = N * w;//Why times jacob again?

    /* eq. (29) */
    auto const  gamma_t_int = jacob * omega_z_v * dv_do_dz_v * zo_zt_v * hmf_v * mor_v * w * lc_lt_v;

    /* eq. (30) */
    auto gamma_t_cen = [this, lnM, zt, A, lo_lt_v, lc_jacob](double radius) {
        return fcen_ * del_sig_cen(radius, lnM, zt) * exp(A * T_cen(radius, lnM)) * lo_lt_v / (lc_jacob * 6.28318530718);
    };

    /* eq. (31) */
    auto gamma_t_mis = [this, lnM, zt, A, R, lo, lc, theta, lc_lt_v](double radius) {
        return (1.0 - fcen_) * roffset(R) * lo_lc(lo, lc, R)
            * del_sig_cen(std::sqrt(radius*radius + R*R + 2*R*radius * std::cos(theta)), lnM, zt)
            * exp(A * T_cen(radius, lnM)) * lc_lt_v / (6.28318530718);
    };

    /* eq. (28) */
    auto const  gamma_t = y3_cluster::transform(r, 
		    [m_shear, Nw, sig_crit, gamma_t_int, gamma_t_cen, gamma_t_mis]
                    (double radius) {
                        /* Nw intentionally left out - returned in return_arr to be used further on */
                        return (1.0 + m_shear) / sig_crit
                                * gamma_t_int * (gamma_t_cen(radius) + gamma_t_mis(radius));
                    });

    std::array<double, NRADII+2> return_arr;
    std::copy_n( gamma_t.begin(), gamma_t.size(),  return_arr.begin() );
    return_arr[NRADII] = N;
    return_arr[NRADII+1] = Nw;


    // this is for y1 likelihood
    //double const N = omega_z_v * dv_do_dz_v * zo_zt_v * hmf_v * mor_v * lc_lt_v ;//*(fcen_+ (1.0-fcen_)*roffset(R)*lo_lc(lo, lc, R));
    //double const Nw = jacob * N * w;
    //double const gamma_t_int = omega_z_v * dv_do_dz_v * zo_zt_v * hmf_v * mor_v * w * lc_lt_v;
    //double const gamma_t_cen =del_sig_cen(r, lnM);//fcen_ * del_sig_cen(r, lnM) * exp(A * T_cen(r, lnM));
    //double const gamma_t_mis = 0.;//(1.0 - fcen_) * roffset(R)*lo_lc(lo, lc, R) * del_sig_mis(r, lnM, R) *  exp(A * T_cen(r, lnM));
    //gamma_t=(1.0 + m_shear)/sig_crit_inv * gamma_t_int * (gamma_t_cen + gamma_t_mis)
    //double const jacob=lnM_ir_.jacobian() * lo_ir_.jacobian() * lt_ir_.jacobian() * lc_ir_.jacobian() * zo_ir_.jacobian() * zt_ir_.jacobian() * R_ir_.jacobian() * A_ir_.jacobian();
    //return {N*jacob, Nw*jacob, gamma_t*jacob}
    return return_arr;
  }
};

template <typename MOR,
          typename LO_LC,
          typename LC_LT,
          typename ZO_ZT,
          typename ROFFSET,
          typename T_CEN,
          typename T_MIS,
          typename A_CEN,
          typename A_MIS,
          typename HMF,
          typename DEL_SIG_CEN,
          /*typename DEL_SIG_MIS,*/
          typename DV_DO_DZ,
          typename OMEGA_Z>
Gamma_T_Integrand<MOR,
                  LO_LC,
                  LC_LT,
                  ZO_ZT,
                  ROFFSET,
                  T_CEN,
                  T_MIS,
                  A_CEN,
                  A_MIS,
                  HMF,
                  DEL_SIG_CEN,
                  /*DEL_SIG_MIS,*/
                  DV_DO_DZ,
                  OMEGA_Z>
make_gamma_t_integrand(double fcen,
                       double msci,
                       MOR mor,
                       LO_LC lo_lc,
                       LC_LT lc_lt,
                       ZO_ZT zo_zt,
                       ROFFSET roffset,
                       T_CEN t_cen,
                       T_MIS t_mis,
                       A_CEN a_cen,
                       A_MIS a_mis,
                       HMF hmf,
                       DEL_SIG_CEN del_sig_cen,
                       /*DEL_SIG_MIS del_sig_mis,*/
                       DV_DO_DZ dv_do_dz,
                       OMEGA_Z omega_z,
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

   std::array<double, 10> rarray; // can I pass a vector here?
   for ( std::size_t i = 0; i < 10; i++ ) {rarray[ i ] = 0.1*(i+0.1);}
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
          /*del_sig_mis,*/
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
