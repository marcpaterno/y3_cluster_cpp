#include <cubacpp/gsl.hh>
#include <iostream>

#include "catch2/catch.hpp"
#include "primitives.hh"
#include "sin_cos_polynomial_integral.hh"

using y3_cluster::sinusoid_polynomial_integral,
      y3_cluster::sinusoid_polynomial_integrals,
      y3_cluster::integer_pow;

TEST_CASE("Test analytic sin/cos polynomial integrals")
{
    cubacpp::QAG qag;

    std::vector<std::pair<double, double>> ranges{{{-7, -3}, {-3, -2}, {-2, -1},
                                                   {-0.3, -0.1}, {0.1, 2},
                                                   {1, 2}, {2, 3}, {3, 7}}};
    const int minpow = -12, maxpow = 12;
    for (const auto [minx, maxx] : ranges) {
        const auto [sin_results, cos_results] = sinusoid_polynomial_integrals(minpow, maxpow, minx, maxx);
        for (auto pow = minpow; pow < maxpow; pow++) {
            const auto sin_integral = qag.with_range(minx, maxx)
                                         .integrate([&](double x) { return integer_pow(x, pow) * sin(x); },
                                                    1e-5, 1e-18);
            const auto cos_integral = qag.with_range(minx, maxx)
                                         .integrate([&](double x) { return integer_pow(x, pow) * cos(x); },
                                                    1e-5, 1e-18);

            REQUIRE(sin_integral.status == 0);
            REQUIRE(cos_integral.status == 0);

            const auto [sin_result, cos_result] = sinusoid_polynomial_integral(pow, minx, maxx);

            CHECK(sin_integral.value == Approx(sin_result).epsilon(1e-5));
            CHECK(cos_integral.value == Approx(cos_result).epsilon(1e-5));

            CHECK(sin_integral.value == Approx(sin_results[pow - minpow]).epsilon(1e-5));
            CHECK(cos_integral.value == Approx(cos_results[pow - minpow]).epsilon(1e-5));
        }
    }

    std::vector<std::pair<double, double>> small_ranges{{{-2.6, -2.5}, {-2.5, -2}, {-2, -0.01},
                                                         {-1, -0.1}, {-0.3, -0.01},
                                                         {0.001, 0.01}, {0.05, 2.3}, {0.1, 0.3}, {1, 2}}};
    for (const auto [minx, maxx] : small_ranges) {
        for (auto pow = minpow; pow < maxpow; pow++) {
            const auto sin_integral = qag.with_range(minx, maxx)
                                         .integrate([&](double x) { return integer_pow(x, pow) * sin(x); },
                                                    1e-5, 1e-18);
            const auto cos_integral = qag.with_range(minx, maxx)
                                         .integrate([&](double x) { return integer_pow(x, pow) * cos(x); },
                                                    1e-5, 1e-18);
            REQUIRE(sin_integral.status == 0);
            REQUIRE(cos_integral.status == 0);

            const auto [sin_result, cos_result] = y3_cluster::sinusoid_polynomial_taylor_approx(pow, minx, maxx);

            CHECK(sin_integral.value == Approx(sin_result).epsilon(1e-4));
            CHECK(cos_integral.value == Approx(cos_result).epsilon(1e-4));
        }
    }
}
