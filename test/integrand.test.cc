#include "catch2/catch.hpp"
#include "cubacpp/integrand.hh"

TEST_CASE("Scalar 1D integrand can be called", "[integrand]")
{
  auto f1 = [](double x) { return 2. * x; };
  auto i1 = cubacpp::detail::integrand<decltype(f1)>;
  double res = 0.0;
  int const ndim = 1;
  double const x = 3.5;
  int const ncomp = 1;
  i1(&ndim, &x, &ncomp, &res, &f1);
  REQUIRE(res == 7.0);
}

TEST_CASE("Scalar 2D integrand can be called", "[integrand]")
{
  auto f2 = [](double x, double y) { return x + y; };
  auto i2 = cubacpp::detail::integrand<decltype(f2)>;
  double res = 0.0;
  int const ndim = 2;
  double const x[] = {2.0, -3.0};
  int const ncomp = 1;
  i2(&ndim, x, &ncomp, &res, &f2);
  REQUIRE(res == -1.0);
}

TEST_CASE("2D vector 1D integrand can be called", "[integrand]")
{
  auto f = [](double x) { return std::array<double, 2>{{-x, x}}; };
  auto igrand = cubacpp::detail::integrand<decltype(f)>;
  double res[] = {0.0, 0.0};
  REQUIRE(sizeof(res) == 2 * sizeof(double));
  double const x = 5.0;
  int const ndim = 1;
  int const ncomp = 2;
  igrand(&ndim, &x, &ncomp, res, &f);
  REQUIRE(res[0] == -5.0);
  REQUIRE(res[1] == 5.0);
}

TEST_CASE("2D vector 2D integrand can be called", "[integrand]")
{
  auto f = [](double x, double y, double z) {
    return std::array<double, 2>{{x + y, y + z}};
  };
  auto igrand = cubacpp::detail::integrand<decltype(f)>;
  double res[] = {0.0, 0.0};
  double const x[] = {1.0, 2.0, 3.0};
  int const ndim = sizeof(x) / sizeof(double);
  int const ncomp = sizeof(res) / sizeof(double);
  igrand(&ndim, x, &ncomp, res, &f);
  REQUIRE(res[0] == 3.0);
  REQUIRE(res[1] == 5.0);
}
