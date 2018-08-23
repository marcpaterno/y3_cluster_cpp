#include <cubacpp/cubacpp.hh>
#include <iostream>

#include "catch2/catch.hpp"
#include "primitives.hh"
#include "sin_cos_polynomial_integral.hh"

using y3_cluster::sinusoid_polynomial_integral,
      y3_cluster::integer_pow;

TEST_CASE("Test analytic sin/cos polynomial integrals")
{
    cubacpp::QAG qag;

    std::vector<std::pair<double, double>> ranges{{{-7, -3}, {-3, -2}, {-2, -1}, {1, 2}, {2, 3}, {3, 7}}};
    for (auto pow = -10; pow < 12; pow++) {
        for (const auto [minx, maxx] : ranges) {
            const auto cos_integral = qag.with_range(minx, maxx)
                                         .integrate([&](double x) { return integer_pow(x, pow) * cos(x); },
                                                    1e-5, 1e-18);
            const auto sin_integral = qag.with_range(minx, maxx)
                                         .integrate([&](double x) { return integer_pow(x, pow) * sin(x); },
                                                    1e-5, 1e-18);

            REQUIRE(cos_integral.status == 0);
            REQUIRE(sin_integral.status == 0);

            const auto [sin_result, cos_result] = sinusoid_polynomial_integral(pow, minx, maxx);

            CHECK(sin_integral.value == Approx(sin_result).epsilon(1e-5));
            CHECK(cos_integral.value == Approx(cos_result).epsilon(1e-5));
        }
    }
}
