#include <cubacpp/cubacpp.hh>
#include <iostream>

#include "catch2/catch.hpp"
#include "primitives.hh"
#include "bessel_polynomial_integral.hh"

using y3_cluster::bessel_polynomial_integral,
      y3_cluster::bessel_polynomial_integrals,
      y3_cluster::sinusoid_polynomial_integral,
      y3_cluster::integer_pow;

TEST_CASE("Test analytic sin/cos polynomial integrals")
{
  cubacpp::QAG qag;

  const auto test_integral = qag.with_range(1, 2)
                                .integrate([](double x) {
                                        return x * gsl_sf_bessel_jl(3, x);
                                        }, 1e-5, 1e-18);

  REQUIRE(test_integral.status == 0);

  const auto test_analytic = [](double x) {
    return -1 * (3 * sinusoid_polynomial_integral(-3, 0.5, x).first
               + (3 / x) * gsl_sf_bessel_jl(0, x)
               + 3 * gsl_sf_bessel_jl(1, x)
               + x * gsl_sf_bessel_jl(2, x));
  };

  CHECK(test_integral.value == Approx(test_analytic(2) - test_analytic(1)));

  std::vector<std::pair<double, double>> ranges{{{1, 2}, {2, 3}, {3, 7}, {10, 20}}};
  std::vector<double> ks{{1.0, 2.0, 5.0, 10.0}};
  const double maxl = 10;
  for (auto pow = 0; pow < 10; pow++) {
    for (const auto [minx, maxx] : ranges) {
      for (const auto k : ks) {
        const auto list = bessel_polynomial_integrals(pow, maxl, k, minx, maxx);
        for (auto l = 0u; l < maxl; l++) {
          const auto bessel_integral = qag.with_range(minx, maxx)
                                          .integrate([&](double x) {
                                                     return integer_pow(x, pow) * gsl_sf_bessel_jl(l, k * x);
                                                  },
                                                  1e-5, 1e-18);

          REQUIRE(bessel_integral.status == 0);

          const auto result = bessel_polynomial_integral(pow, l, k, minx, maxx);

          if (!(bessel_integral.value == Approx(result).epsilon(1e-3).margin(1e-10)))
              std::cout << "x^" << pow << " * j_" << l << "(" << k << " x) on ["
                        << minx << ", " << maxx << "]\n";
          CHECK(bessel_integral.value == Approx(result).epsilon(1e-3).margin(1e-10));
          CHECK(bessel_integral.value == Approx(list[l]).epsilon(1e-3).margin(1e-10));
        }
      }
    }
  }
}
