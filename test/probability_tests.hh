#include "catch2/catch.hpp"
#include "test/lc_lt_t.hh"
#include "test/lo_lc_t.hh"
#include "test/roffset_t.hh"
#include "cubacpp/cubacpp.hh"
#include "integration_range.hh"

#include <iostream>

namespace y3_cluster {

    template<typename Integrator>
    void
    test_integrate_lc_lt(Integrator I,
                         IntegrationRange lc_ir,
                         IntegrationRange lt_ir,
                         IntegrationRange zt_ir,
                         std::size_t lt_width,
                         std::size_t zt_width,
                         const double epsrel=1.0e-4,
                         const double epsabs=1.0e-12,
                         bool print=false,
                         bool test=true)
    {
        if (print)
            std::cout << "lt,zt,lc_lt_integrated,status,error,prob\n";

        LC_LT_t lc_lt;
        for (auto i = 0u; i < lt_width; i++) {
            for (auto j = 0u; j < zt_width; j++) {
                const double lt = lt_ir.transform((i + 1) / ((double) lt_width + 1));
                const double zt = zt_ir.transform((j + 1) / ((double) zt_width + 1));
                const auto res = I.integrate([lt, zt, lc_ir, lc_lt](double scaled_lc) {
                          const double lc = lc_ir.transform(scaled_lc);
                          return lc_ir.jacobian() * lc_lt(lc, lt, zt);
                        },
                        epsrel, epsabs);

                if (print)
                    std::cout << lt << ", "
                              << zt << ", "
                              << res.value << ", "
                              << res.status << ", "
                              << res.error << ", "
                              << res.prob << '\n';

                if (test) {
                    // In reality - since we are integrating [1, 200] not [0, +inf] - it
                    // won't be exactly 1.0. So, arbitrary wiggle room 1e-2
                    // TODO: should this epsilon be changed?
                    CHECK(res.status == 0);
                    CHECK(res.value == Approx(1.0).epsilon(1e-2).margin(1e-2));
                }
            }
        }
    }

    template<typename Integrator>
    void
    test_integrate_lo_lc(Integrator I,
                         IntegrationRange lo_ir,
                         IntegrationRange lc_ir,
                         IntegrationRange R_ir,
                         std::size_t lc_width,
                         std::size_t R_width,
                         const double epsrel=1.0e-4,
                         const double epsabs=1.0e-12,
                         bool print=false,
                         bool test=true)
    {
        if (print)
            std::cout << "lc,R,lc_lt_integrated,status,error,prob\n";

        // Values from trivial_gamma_t
        LO_LC_t lo_lc{1.66, 0.26, 1.43, 1.0};

        for (auto i = 0u; i < lc_width; i++) {
            for (auto j = 0u; j < R_width; j++) {
                const double lc = lc_ir.transform((i + 1) / ((double) lc_width + 1));
                const double R = R_ir.transform((j + 1) / ((double) R_width + 1));
                const auto res = I.integrate([lc, R, lo_ir, lo_lc](double scaled_lo) {
                          const double lo = lo_ir.transform(scaled_lo);
                          return lo_ir.jacobian() * lo_lc(lo, lc, R);
                        },
                        epsrel, epsabs);

                if (print)
                    std::cout << lc << ", "
                              << R << ", "
                              << res.value << ", "
                              << res.status << ", "
                              << res.error << ", "
                              << res.prob << '\n';

                if (test) {
                    CHECK(res.status == 0);
                    CHECK(res.value == Approx(1.0).epsilon(1e-2).margin(5e-2));
                }
            }
        }
    }

    template<typename Integrator>
    void
    test_integrate_roffset(Integrator I,
                           IntegrationRange R_ir,
                           double tau,
                           const double epsrel=1.0e-4,
                           const double epsabs=1.0e-12,
                           bool print=false,
                           bool test=true)
    {
        if (print)
            std::cout << "roffset_integrand,status,error,prob\n";

        // Values from trivial_gamma_t
        ROFFSET_t roffset(tau);

        const auto res = I.integrate([R_ir, roffset](double scaled_R) {
                  const double R = R_ir.transform(scaled_R);
                  return R_ir.jacobian() * roffset(R);
                },
                epsrel, epsabs);

        if (print)
            std::cout << res.value << ", "
                      << res.status << ", "
                      << res.error << ", "
                      << res.prob << '\n';

        if (test) {
            CHECK(res.status == 0);
            CHECK(res.value == Approx(1.0).epsilon(1e-2).margin(5e-2));
        }
    }

}
