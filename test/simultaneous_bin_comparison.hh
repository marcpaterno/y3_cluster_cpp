#include "catch2/catch.hpp"
#include "gamma_t.hh"

#include <algorithm>


namespace y3_cluster {

    template<typename Integrator, typename MODELS,
             std::size_t NRADII, std::size_t NRICHNESS, std::size_t NREDSHIFT>
    void
    test_simultaneous_bins(Integrator I,
                           Gamma_T_Integrand<MODELS, NRADII, NRICHNESS, NREDSHIFT> gti,
                           const double epsrel=1e-4,
                           const double epsabs=1e-12)
    {
        auto res_cen = I.integrate([&gti](double a, double b, double c,
                                          double d, double e) {
                                        return gti.centered(a, b, c, d, e);
                                    },
                                    epsrel, epsabs);

        CHECK(res_cen.status == 0);

        auto res_mis = I.integrate([&gti](double a, double b, double c,
                                          double d, double e, double f,
                                          double g, double h) {
                                        return gti.miscentered(a, b, c, d, e, f, g, h);
                                    },
                                    epsrel, epsabs);

        CHECK(res_mis.status == 0);

        auto bins_cen = make_gamma_t_integrated_bins(gti, res_cen);
        auto bins_mis = make_gamma_t_integrated_bins(gti, res_mis);

        for (auto i = 0u; i < NRICHNESS * NREDSHIFT; i++) {
            const auto& bin_cen = bins_cen[i];
            const auto& bin_mis = bins_mis[i];

            std::array<y3_cluster::IntegrationRange, 1> new_lir = {bin_cen.lo_ir};
            std::array<y3_cluster::IntegrationRange, 1> new_zir = {bin_cen.zo_ir};
            auto gti_new = gti.with_bins(new_lir, new_zir);
            auto checker = [epsabs, epsrel](const double a, const double b) {
                CHECK(a == Approx(b).epsilon(2 * epsrel).margin(2 * epsabs));
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
        }
    }

}
