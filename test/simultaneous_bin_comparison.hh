#include "catch2/catch.hpp"
#include "gamma_t.hh"

#include <algorithm>
#include <iostream>

namespace y3_cluster {

    template<typename Integrator, typename MODELS,
             std::size_t NRADII, std::size_t NRICHNESS, std::size_t NREDSHIFT>
    void
    test_simultaneous_bins(Integrator I,
                           Gamma_T_Integrand<MODELS, NRADII, NRICHNESS, NREDSHIFT> gti,
                           const double epsrel=1e-3,
                           const double epsabs=1e-12,
                           bool print = false,
                           bool test = true)
    {
        auto res_cen = I.integrate([&gti](double a, double b, double c,
                                          double d, double e) {
                                        return gti.centered(a, b, c, d, e);
                                    },
                                    epsrel, epsabs);

        if (test)
            CHECK(res_cen.status == 0);

        auto res_mis = I.integrate([&gti](double a, double b, double c,
                                          double d, double e, double f,
                                          double g, double h) {
                                        return gti.miscentered(a, b, c, d, e, f, g, h);
                                    },
                                    epsrel, epsabs);

        if (test)
            CHECK(res_mis.status == 0);

        auto bins_cen = make_gamma_t_integrated_bins(gti, res_cen);
        auto bins_mis = make_gamma_t_integrated_bins(gti, res_mis);

        // CSV Output Header
        if (print)
            std::cout << "lo_min,lo_max,zo_min,zo_max,"
                      << "cen_or_mis,"
                      << "N_single,N_binned,N_single_status,N_binned_status,"
                      << "N_single_error,N_binned_error,N_single_prob,N_binned_prob"
                      << '\n';

        for (auto i = 0u; i < NRICHNESS * NREDSHIFT; i++) {
            const auto& bin_cen = bins_cen[i];
            const auto& bin_mis = bins_mis[i];

            std::array<y3_cluster::IntegrationRange, 1> new_lir = {bin_cen.lo_ir};
            std::array<y3_cluster::IntegrationRange, 1> new_zir = {bin_cen.zo_ir};
            auto gti_new = gti.with_bins(new_lir, new_zir);
            auto checker = [test, epsabs, epsrel](const double a, const double b) {
                const double rel_diff = 2 * epsrel * std::max(a, b);
                if (test)
                    CHECK(a == Approx(b).margin(std::max(rel_diff, 2 * epsabs)));
            };

            auto print_line = [print, bin_cen](const char *cen_or_mis,
                                               const auto bin, const auto bin_single,
                                               const auto res, const auto res_single) {
                if (print)
                    std::cout << bin_cen.lo_ir.transform(0.0) << ", "
                              << bin_cen.lo_ir.transform(1.0) << ", "
                              << bin_cen.zo_ir.transform(0.0) << ", "
                              << bin_cen.zo_ir.transform(1.0) << ", "
                              << cen_or_mis << ","
                              << bin_single.N << ", "
                              << bin.N << ", "
                              << res_single.status << ", "
                              << res.status << ", "
                              << bin_single.N_error << ", "
                              << bin.N_error << ", "
                              << bin_single.N_prob << ", "
                              << bin.N_prob << '\n';
            };

            // First, check that the centered term gives comparable results
            auto res_cen_single = I.integrate([&gti_new](double a, double b, double c,
                                                         double d, double e) {
                                                   return gti_new.centered(a, b, c, d, e);
                                               },
                                               epsrel, epsabs);

            CHECK(res_cen_single.status == 0);

            const auto bin_cen_single = make_gamma_t_integrated_bins(gti_new, res_cen_single)[0];
            for(auto j = 0u; j < NRADII; j++)
                checker(bin_cen_single.gamma_ts[i], bin_cen.gamma_ts[i]);

            checker(bin_cen_single.N, bin_cen.N);
            checker(bin_cen_single.Nw, bin_cen.Nw);

            print_line("cen", bin_cen, bin_cen_single, res_cen, res_cen_single);

            // Now test the miscentered term
            auto res_mis_single = I.integrate([&gti_new](double a, double b, double c,
                                                         double d, double e, double f,
                                                         double g, double h) {
                                                       return gti_new.miscentered(a, b, c, d, e, f, g, h);
                                                   },
                                                   epsrel, epsabs);

            CHECK(res_mis_single.status == 0);

            const auto bin_mis_single = make_gamma_t_integrated_bins(gti_new, res_mis_single)[0];
            for(auto j = 0u; j < NRADII; j++)
                checker(bin_mis_single.gamma_ts[i], bin_mis.gamma_ts[i]);

            checker(bin_mis_single.N, bin_mis.N);
            checker(bin_mis_single.Nw, bin_mis.Nw);

            print_line("mis", bin_mis, bin_mis_single, res_mis, res_mis_single);
        }
    }

}
