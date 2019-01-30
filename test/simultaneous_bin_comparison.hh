#include "catch2/catch.hpp"
#include "gamma_t.hh"

#include <algorithm>
#include <iostream>

namespace y3_cluster {

    template<typename Integrator, typename MODELS>
    void
    test_simultaneous_bins(Integrator I,
                           Gamma_T_Integrand<MODELS>  gti,
                           const double epsrel=1e-3,
                           const double epsabs=1e-12,
                           bool print = false,
                           bool test = true)
    {
        auto bins_cen = gti.integrate_centered(I, epsrel, epsabs);
        if (test)
            CHECK(bins_cen.status == 0);

        auto bins_mis = gti.integrate_miscentered(I, epsrel, epsabs);
        if (test)
            CHECK(bins_mis.status == 0);

        // CSV Output Header
        if (print)
            std::cout << "lo_min,lo_max,zo_min,zo_max,"
                      << "cen_or_mis,"
                      << "N_single,N_binned,N_single_status,N_binned_status,"
                      << "N_single_error,N_binned_error,N_single_prob,N_binned_prob"
                      << '\n';

        auto const top = gti.lo_ir_.size () * gti.zo_ir_.size ();

        for (auto i = 0u; i < top; i++) {
            const auto& bin_cen = bins_cen[i];
            const auto& bin_mis = bins_mis[i];

            std::vector<y3_cluster::IntegrationRange> new_lir = {bin_cen.lo_ir};
            std::vector<y3_cluster::IntegrationRange> new_zir = {bin_cen.zo_ir};
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
            auto bin_cen_single = gti_new.integrate_centered(I, epsrel, epsabs);
            CHECK(bin_cen_single.status == 0);

            for (auto j = 0u; j < bin_cen_single[0].gamma_ts.size(); j++)
                for (auto k = 0u; k < bin_cen_single[0].gamma_ts[0].size(); k++)
                    checker(bin_cen_single[0].gamma_ts[j][k], bin_cen.gamma_ts[j][k]);

            checker(bin_cen_single[0].N, bin_cen.N);
            checker(bin_cen_single[0].Nw, bin_cen.Nw);

            print_line("cen", bin_cen, bin_cen_single[0], bins_cen, bin_cen_single);

            // Now test the miscentered term
            auto bin_mis_single = gti_new.integrate_miscentered(I, epsrel, epsabs);
            CHECK(bin_mis_single.status == 0);

            for (auto j = 0u; j < bin_mis_single[0].gamma_ts.size(); j++)
                for (auto k = 0u; k < bin_mis_single[0].gamma_ts[0].size(); k++)
                    checker(bin_mis_single[0].gamma_ts[j][k], bin_mis.gamma_ts[j][k]);

            checker(bin_mis_single[0].N, bin_mis.N);
            checker(bin_mis_single[0].Nw, bin_mis.Nw);

            print_line("mis", bin_mis, bin_mis_single[0], bins_mis, bin_mis_single);
        }
    }

}
