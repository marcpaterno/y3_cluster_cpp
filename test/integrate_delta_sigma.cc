#include <iostream>
#include <cubacpp/cubacpp.hh>
#include <integration_range.hh>

#include <test/lc_lt_t.hh>
#include <test/lo_lc_t.hh>
#include <test/param_space_explorer.hh>
#include <test/t_cen_t.hh>
#include <test/roffset_t.hh>
#include <test/del_sig_cen_y1.hh>

struct delta_sigma_parameters {
    double lo, lt, zt, R, lnM;

    delta_sigma_parameters(double lo, double lt, double zt, double R, double lnM)
        : lo(lo)
        , lt(lt)
        , zt(zt)
        , R(R)
        , lnM(lnM)
        {}
};

using y3_cluster::LO_LC_t,
      y3_cluster::LC_LT_t,
      y3_cluster::ROFFSET_t,
      y3_cluster::T_CEN_t,
      y3_cluster::DEL_SIG_CEN_y1,
      y3_cluster::IntegrationRange;

template<typename Integrator>
cubacpp::integration_results<1>
integrate_del_sig_cen(Integrator I,
                      delta_sigma_parameters& params,
                      LC_LT_t lc_lt,
                      T_CEN_t tcen,
                      DEL_SIG_CEN_y1 dsc,
                      IntegrationRange A_range,
                      double epsrel = 1e-3,
                      double epsabs = 1e-12)
{
    const auto& lo = params.lo,
                lt = params.lt,
                zt = params.zt,
                R = params.R,
                lnM = params.lnM;

    const double lo_lt_v = lc_lt(lo, lt, zt);
    return I.integrate([&](double scaled_A) {
                const double A = A_range.transform(scaled_A);
                return lo_lt_v
                       * exp(A * tcen(R, lnM))
                       * dsc(R, lnM, zt);
            },
            epsrel,
            epsabs);
}


template<typename Integrator>
cubacpp::integration_results<1>
integrate_del_sig_mis(Integrator I,
                      delta_sigma_parameters& params,
                      LC_LT_t lc_lt,
                      LO_LC_t lo_lc,
                      ROFFSET_t roffset,
                      T_CEN_t tcen,
                      DEL_SIG_CEN_y1 dsc,
                      IntegrationRange R_mis_range,
                      IntegrationRange theta_range,
                      IntegrationRange lc_range,
                      IntegrationRange A_range,
                      double epsrel = 1e-3,
                      double epsabs = 1e-12)
{
    const auto& lo = params.lo,
                lt = params.lt,
                zt = params.zt,
                R = params.R,
                lnM = params.lnM;

    return I.integrate([&](double scaled_R_mis,
                           double scaled_theta,
                           double scaled_lc,
                           double scaled_A) {
                const double R_mis = R_mis_range.transform(scaled_R_mis),
                             theta = theta_range.transform(scaled_theta),
                             lc = lc_range.transform(scaled_lc),
                             A = A_range.transform(scaled_A),
                             lc_lt_v = lc_lt(lc, lt, zt),
                             R_adj = sqrt(R*R + R_mis*R_mis + 2*R*R_mis * cos(theta));
                return lc_lt_v
                       * lo_lc(lo, lc, R_mis)
                       * roffset(R_mis)
                       * exp(A * tcen(R_adj, lnM))
                       * dsc(R_adj, lnM, zt);
            },
            epsrel,
            epsabs);
}

int
main()
{
    std::cout << "lo,lt,zt,R,lnM"
              << ",del_sig_cen,cen_status,cen_error,cen_prob"
              << ",del_sig_mis,mis_status,mis_error,mis_prob"
              << '\n';

    delta_sigma_parameters dsp(25, 25,
                               0.2,
                               0.5,
                               std::log(1e15));

    // Parameter exploration ranges
    IntegrationRange lo_ir{20, 28},
                     lt_ir{2, 50},
                     zt_ir{0.1, 0.3},
                     R_ir{0., 1.0},
                     lnM_ir{std::log(5e11), std::log(1e16)};

    // Proper integration ranges
    IntegrationRange A_ir{-0.01, 0.01},
                     R_i{0.1, 1},
                     theta_ir{0.0, 6.28318530718},
                     lc_ir{2, 50};

    ParamSpaceExplorer<5> pse{{&dsp.lo, &dsp.lt, &dsp.zt, &dsp.R, &dsp.lnM},
                              {lo_ir, lt_ir, zt_ir, R_ir, lnM_ir},
                              {3, 3, 3, 3, 3}};

    LO_LC_t lo_lc{1.66, 0.26, 1.43, 1.0};
    LC_LT_t lc_lt;
    ROFFSET_t roffset{0.2};
    T_CEN_t tcen;
    DEL_SIG_CEN_y1 dsc;

    cubacpp::Vegas c;
    c.maxeval = 999999999;

    do {
        auto res_cen = integrate_del_sig_cen(c,
                                             dsp,
                                             lc_lt,
                                             tcen,
                                             dsc,
                                             A_ir);
        auto res_mis = integrate_del_sig_mis(c,
                                             dsp,
                                             lc_lt,
                                             lo_lc,
                                             roffset,
                                             tcen,
                                             dsc,
                                             R_ir,
                                             theta_ir,
                                             lc_ir,
                                             A_ir);
        std::cout << dsp.lo << "," << dsp.lt << "," << dsp.zt << "," << dsp.R << "," << dsp.lnM << ","
                  << res_cen.value << "," << res_cen.status << "," << res_cen.error << "," << res_cen.prob << ","
                  << res_mis.value << "," << res_mis.status << "," << res_mis.error << "," << res_mis.prob
                  << '\n';
    } while (pse.next());

    return 0;
}
