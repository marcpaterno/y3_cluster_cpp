#include <cubacpp/cubacpp.hh>
#include <iostream>

#include "catch2/catch.hpp"
#include "sin_cos_polynomial_integrals.hh"

using y3_cluster::sinusoid_polynomial_integral;

constexpr double
integer_pow(double n, int pow)
{
    if (pow == 0)
        return 1;
    if (pow < 0)
        return integer_pow(n, pow + 1) / n;
    return integer_pow(n, pow - 1) * n;
}

TEST_CASE("Test analytic sin/cos polynomial integrals")
{
    cubacpp::QAG qag;

    CHECK(integer_pow(2, -1) == Approx(0.5));
    CHECK(integer_pow(3, -2) == Approx(1.0/3.0/3.0));

    std::vector<std::pair<double, double>> ranges{{{1, 2}, {2, 3}, {3, 7}}};
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
