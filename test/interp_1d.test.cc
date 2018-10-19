#include "catch2/catch.hpp"
#include "interp_1d.hh"
#include "transform.hh"

#include <cmath>
#include <cstddef>
#include <stdexcept>

using y3_cluster::Interp1D;

TEST_CASE("Interp1D exact at knots", "[interpolation][1d]")
{
  std::vector<double> const xs = {1., 2., 3., 4., 5., 6, 7., 8., 9};
  std::vector<double> const ys
    = y3_cluster::transform<double,double>
              ( xs,    [](double x) { return 2 * x * (3 - x) * std::cos(x); } );
  Interp1D f{xs, ys};
  for (std::size_t i = 0; i != xs.size (); ++i) {
    CHECK(ys[i] == f(xs[i]));
  }
}

TEST_CASE("Interp1D on quadratic")
{
  std::vector<double> const xs = {1., 2., 3.};
  auto fcn = [](double x) { return x * x; };
  std::vector<double> const ys = y3_cluster::transform<double,double>(xs, fcn);
  Interp1D f{xs, ys};
  // Answer produced using Mathematica 11.3.0.0,
  // using  Interpolation[{1., 4., 9.}, InterpolationOrder -> 1]
  CHECK(f(1.41421) == Approx(2.24263).epsilon(1e-4));
}

TEST_CASE("Interp1D throws on extrapolation")
{
  std::vector<double> const xs = {-5., 5.};
  auto fcn = [](double x) { return -4 * x; };
  std::vector<double> const ys = y3_cluster::transform<double,double>(xs, fcn);
  Interp1D f{xs, ys};
  CHECK_THROWS_AS(f(-10), std::domain_error);
  CHECK_NOTHROW(f(-5.0));
  CHECK_NOTHROW(f(5.0));
  CHECK_THROWS_AS(f(10), std::domain_error);
}
