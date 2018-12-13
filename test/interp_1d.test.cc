#include "catch2/catch.hpp"
#include "interp_1d.hh"
#include "transform.hh"

#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <stdexcept>

using y3_cluster::Interp1D;

TEST_CASE("Interp1D exact at knots", "[interpolation][1d]")
{
  std::vector<double> const xs = {1., 2., 3., 4., 5., 6, 7., 8., 9};
  std::vector<double> const ys
    = y3_cluster::transform<double,double>
              ( xs,    [](double x) { return 2 * x * (3 - x) * std::cos(x); } );
  Interp1D f{xs, ys}; for (std::size_t i = 0; i != xs.size (); ++i) {
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

TEST_CASE("Validate")
{
  std::array const xs { 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};
  auto fcn =[](double x){ return x*x; };
  std::array const ys = y3_cluster::transform(xs, fcn);
  Interp1D f{xs, ys};
  using std::sqrt;
  std::array const xtest { sqrt(1.1), sqrt(1.2), sqrt(1.3), sqrt(1.4),
                           sqrt(2.5), sqrt(2.6), sqrt(2.7), sqrt(2.8),
                           sqrt(3.9) };
  std::ofstream out { "../data/interp_1d.out" };
  out << "x\tytrue\tytest\n";
  for (auto x : xtest)
  {
    out << x << '\t' << fcn(x) << '\t' << f(x) << '\n';
  }
}
