#include "catch2/catch.hpp"
#include "cubacpp/cubacpp.hh"
#include "utils/integration_range.hh"
#include <array>
#include <cmath>
#include <iostream>
#include <tuple>

// MyFunc is an example of a user-defined function to be integrated. It to be
// written as either a class or a struct. It must have const member function
// operator(), the function call operator, which takes one or more doubles (or
// types that can be converted to doubles).

double constexpr pi = 0x1.921fb54442d18p+1;
double sf6integral = 4.14199; // From Mathematica 11.2 NIntegrate to 6 figures.
double constexpr epsrel = 1.e-3;
double constexpr epsabs = 1.e-12;

class MyFunc {
public:
  explicit MyFunc(double mul) : multiplier(mul){};

  double
  operator()(double x, double y) const
  {
    return multiplier * x * y * (x + y);
  }

private:
  double multiplier;
};

// Scalar-valued free function of one argument.
inline double
sf1(double x)
{
  auto sinx = std::sin(pi * x);
  return sinx * sinx;
}
double constexpr sf1res = 0.5;

// Scalar-valued free function of two arguments.
inline double
sf2(double x, double y)
{
  return 3. * x * y * (x + y);
}
double constexpr sf2res = 1.0;

// Scalar-value function of two arguments.
// This function represents the same integrand as sf2,
// but has built-in the range of integration of x=a to x=b
// and y=c to y=d.
class sf2_a_b {
public:
  sf2_a_b(double a, double b, double c, double d) : _irx(a, b), _iry(c, d) {}
  sf2_a_b(y3_cluster::IntegrationRange x, y3_cluster::IntegrationRange y)
    : _irx(x), _iry(y)
  {}

  double
  operator()(double x, double y) const
  {
    return _irx.jacobian() * _iry.jacobian() *
           sf2(_irx.transform(x), _iry.transform(y));
  }

private:
  y3_cluster::IntegrationRange _irx;
  y3_cluster::IntegrationRange _iry;
};

// v3f2 is an example vector-valued function of two arguments.
inline std::array<double, 3>
v3f2(double x, double y)
{
  double f1 = x + y;
  double f2 = x * f1;
  double f3 = y * f1;
  return {{f1, f2, f3}};
}
std::array<double, 3> constexpr v3f2res{{1, 7. / 12., 7. / 12.}};

// sf6 is an example vector-valued function of 6 arguments.
inline double
sf6(double u, double v, double w, double x, double y, double z)
{
  return 44 * (std::sin(pi * u * v) / u) * x * y * std::pow(w, y) *
         std::sin(pi * z * z);
}
double sf6res = 4.14199; // From Mathematica 11.2 NIntegrate to 6 figures.

// fracerr returns the absolute value of the fractional error.
double
fracerr(double actual, double expected)
{
  return std::abs((actual - expected) / expected);
}

TEST_CASE("cuhre works", "[integration][cuhre]")
{
  cubacpp::turn_off_cuba_forking();
  cubacpp::Cuhre alg;
  SECTION("sf2")
  {
    auto res = alg.integrate(sf2, epsrel, epsabs);
    CHECK(res.value == Approx(sf2res).epsilon(epsrel));
    CHECK(fracerr(res.value, sf2res) < epsrel);
    CHECK(res.neval < 10000);
    CHECK(res.status == 0);
  }
  SECTION("sf2_a_b")
  {
    sf2_a_b igrand{-4., 7., -11, -3.};
    double sf2_a_b_res = -1276.0;
    auto res = alg.integrate(igrand, epsrel, epsabs);
    CHECK(res.value == Approx(sf2_a_b_res).epsilon(epsrel));
    CHECK(fracerr(res.value, sf2_a_b_res) < epsrel);
    CHECK(res.neval < 10000);
    CHECK(res.status == 0);
  }
  SECTION("v3f2")
  {
    auto res = alg.integrate(v3f2, epsrel, epsabs);
    CHECK(res.value[0] == Approx(1.0).epsilon(epsrel));
    CHECK(res.value[1] == Approx(7. / 12.).epsilon(epsrel));
    CHECK(res.value[2] == Approx(7. / 12.).epsilon(epsrel));
    for (std::size_t i = 0; i != 3; ++i) {
      CHECK(fracerr(res.value[i], v3f2res[i]) < epsrel);
    }
    CHECK(res.status == 0);
  }
  SECTION("myfunc")
  {
    MyFunc ff{3.0};
    auto res = alg.integrate(ff, epsrel, epsabs);
    CHECK(res.value == Approx(1.0).epsilon(epsrel));
    CHECK(fracerr(res.value, 1.) < epsrel);
    CHECK(res.status == 0);
  }
  SECTION("sf6")
  {
    auto res = alg.integrate(sf6, epsrel, epsabs);
    CHECK(res.value == Approx(sf6res).epsilon(epsrel));
    CHECK(fracerr(res.value, sf6res) < epsrel);
    CHECK(res.status == 0);
  }
}

TEST_CASE("vegas works", "[integration][vegas]")
{
  cubacpp::turn_off_cuba_forking();
  cubacpp::Vegas alg;
  SECTION("sf1")
  {
    auto res = alg.integrate(sf1, epsrel, epsabs);
    CHECK(res.value == Approx(0.5).epsilon(epsrel));
    CHECK(fracerr(res.value, sf1res) < epsrel);
    CHECK(res.status == 0);
  }
  SECTION("sf2")
  {
    alg.maxeval = 50 * 1000;
    auto res = alg.integrate(sf2, epsrel, epsabs);
    CHECK(res.value == Approx(1.0).epsilon(epsrel));
    CHECK(fracerr(res.value, sf2res) < epsrel);
    CHECK(res.neval < 15000);
    CHECK(res.status == 0);
  }
  SECTION("v3f2")
  {
    alg.maxeval = 200 * 1000;
    auto res = alg.integrate(v3f2, epsrel, epsabs);
    CHECK(res.value[0] == Approx(1.0).epsilon(epsrel));
    CHECK(res.value[1] == Approx(7. / 12.).epsilon(epsrel));
    CHECK(res.value[2] == Approx(7. / 12.).epsilon(epsrel));
    for (std::size_t i = 0; i != 3; ++i) {
      CHECK(fracerr(res.value[i], v3f2res[i]) < epsrel);
    }
    CHECK(res.status == 0);
  }
  SECTION("myfunc")
  {
    alg.maxeval = 50 * 1000;
    MyFunc ff{3.0};
    auto res = alg.integrate(ff, epsrel, epsabs);
    CHECK(res.value == Approx(1.0).epsilon(epsrel));
    CHECK(fracerr(res.value, 1) < epsrel);
    CHECK(res.status == 0);
  }
  SECTION("sf6")
  {
    alg.maxeval = 100 * 1000;
    auto res = alg.integrate(sf6, epsrel, epsabs);
    CHECK(res.value == Approx(sf6res).epsilon(epsrel));
    CHECK(fracerr(res.value, sf6res) < epsrel);
    CHECK(res.status == 0);
  }
}

TEST_CASE("suave works", "[integration][suave]")
{
  cubacpp::turn_off_cuba_forking();
  cubacpp::Suave alg;
  SECTION("sf1")
  {
    auto res = alg.integrate(sf1, epsrel, epsabs);
    CHECK(res.value == Approx(0.5).epsilon(epsrel));
    CHECK(fracerr(res.value, sf1res) < epsrel);
    CHECK(res.status == 0);
  }
  SECTION("sf2")
  {
    alg.maxeval = 50 * 1000;
    auto res = alg.integrate(sf2, epsrel, epsabs);
    CHECK(res.value == Approx(1.0).epsilon(epsrel));
    CHECK(fracerr(res.value, sf2res) < epsrel);
    CHECK(res.neval < 15000);
    CHECK(res.status == 0);
  }
  SECTION("v3f2")
  {
    alg.maxeval = 50 * 1000;
    auto res = alg.integrate(v3f2, epsrel, epsabs);
    CHECK(res.value[0] == Approx(1.0).epsilon(epsrel));
    CHECK(res.value[1] == Approx(7. / 12.).epsilon(epsrel));
    CHECK(res.value[2] == Approx(7. / 12.).epsilon(epsrel));
    for (std::size_t i = 0; i != 3; ++i) {
      CHECK(fracerr(res.value[i], v3f2res[i]) < epsrel);
    }
    CHECK(res.status == 0);
  }
  SECTION("myfunc")
  {
    alg.maxeval = 50 * 1000;
    MyFunc ff{3.0};
    auto res = alg.integrate(ff, epsrel, epsabs);
    CHECK(res.value == Approx(1.0).epsilon(epsrel));
    CHECK(fracerr(res.value, 1.) < epsrel);
    CHECK(res.status == 0);
  }
  SECTION("sf6")
  {
    alg.maxeval = 100 * 1000;
    auto res = alg.integrate(sf6, epsrel, epsabs);
    CHECK(res.value == Approx(sf6res).epsilon(epsrel));
    CHECK(fracerr(res.value, sf6res) < epsrel);
    CHECK(res.status == 0);
  }
}
