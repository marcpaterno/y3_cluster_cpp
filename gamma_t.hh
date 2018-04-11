#ifndef Y3_CLUSTER_GAMMA_T_HH
#define Y3_CLUSTER_GAMMA_T_HH

#include <cmath>

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
          typename DEL_SIG_MIS,
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
  DEL_SIG_MIS del_sig_mis;
  DV_DO_DZ dv_do_dz;
  OMEGA_Z omega_z;

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
                    DEL_SIG_MIS del_sig_mis,
                    DV_DO_DZ dv_do_dz,
                    OMEGA_Z omega_z)
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
  {}

  // The function call operator -- this is the function to be integrated.
  double
  operator()(double lo,
             double lc,
             double lt,
             double zo,
             double zt,
             double r,
             double R,
             double lnM,
             double A) const
  {
    // We probably should factor out the common subexpressions, rather than
    // relying upon the optimizer to do a perfect job of this for us. This
    // seems to be the intent of the commented-out code below.
    using std::exp;
    auto const hmf_v = hmf(lnM, zt);
    auto const zo_zt_v = zo_zt(zo, zt);
    auto const lc_lt_v = lc_lt(lc, lt, zt);
    auto const mor_v = mor(lt, lnM, zt);
    auto const dv_do_dz_v = dv_do_dz(zt);
    auto const omega_z_v = omega_z(zt);

    // These will eventually be passed by CosmoSIS
    // double omega_zt = 1.0;
    // double dv_do_dz = 1.0;
    double m_shear = 1.0;
    double sig_crit_inv = 1.0;

    // This is the lambda-redshift bin weight that we don't fully understand
    // yet...
    double w = 1.0;

    // The evaluation below follows the convention set in main overleaf document

    double const gamma_t_int =
      omega_z_v * dv_do_dz_v * zo_zt_v * hmf_v * mor_v * w * lc_lt_v;

    double const gamma_t_cen =
      fcen_ * exp(A * T_cen(R, lnM)) * del_sig_cen(r, lnM);

    double const gamma_t_mis = (1.0 - fcen_) * exp(A * T_mis(r, lnM, R)) *
                               del_sig_mis(r, lnM, R) * roffset(R);

    // TODO: Actually calculate Nw. It is itself a multi-dimensional integral
    // for each sampling, so this will take some thought.
    // MFP: If "each sampling" means each sample in the CosmoSIS MCMC, then
    // we can calculate the value of Nw for the current sample, and store
    // that value as a data member in the Gamma_T_Integrand object.
    double const Nw = 1.0;
    double const gamma_t = ((1.0 + m_shear) / (Nw * sig_crit_inv)) *
                           gamma_t_int * (gamma_t_cen + gamma_t_mis);

    return gamma_t;
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
          typename DEL_SIG_MIS,
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
                  DEL_SIG_MIS,
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
                       DEL_SIG_MIS del_sig_mis,
                       DV_DO_DZ dv_do_dz,
                       OMEGA_Z omega_z)
{
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
          del_sig_mis,
          dv_do_dz,
          omega_z};
}

#endif
