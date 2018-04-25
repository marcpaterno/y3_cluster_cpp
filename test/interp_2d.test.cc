#include "catch2/catch.hpp"
#include "test/interp_2d.hh"
#include "test/transform.hh"

#include <cmath>
#include <cstddef>
#include <stdexcept>

using y3_cluster::Interp2D;
template <std::size_t M, std::size_t N>
using matrix = std::array<std::array<double, M>, N>;

TEST_CASE("Interp2D exact at knots", "[interpolation][2d]")
{
  constexpr std::size_t nx = 4;
  constexpr std::size_t ny = 5;
  std::array<double, nx> const xs = {1., 2., 3., 4.};
  std::array<double, ny> const ys = {
    1.,
    2.,
    3.,
    4.,
    5.,
  };
  auto fxy = [](double x, double y) { return 3 * x * y + 2 * x + 4 * y; };
  matrix<4, 5> zs;
  for (std::size_t i = 0; i != nx; ++i) {
    double x = xs[i];
    for (std::size_t j = 0; j != ny; ++j) {
      double y = ys[j];
      zs[i][j] = fxy(x, y);
    }
  }
  Interp2D f{xs, ys, zs};
  for (std::size_t i = 0; i != nx; ++i) {
    double x = xs[i];
    for (std::size_t j = 0; j != ny; ++j) {
      double y = ys[j];
      CHECK(zs[i][j] == f(x, y));
    }
  }
}

// TEST_CASE("Interp2D on quadratic")
// {
//   constexpr const std::size_t numpoints = 3;
//   std::array<double, numpoints> const xs = {1., 2., 3.};
//   auto fcn = [](double x) { return x * x; };
//   std::array<double, numpoints> const ys = y3_cluster::transform(xs, fcn);
//   Interp2D f{xs, ys};
//   // Answer produced using Mathematica 11.3.0.0,
//   // using  Interpolation[{1., 4., 9.}, InterpolationOrder -> 1]
//   CHECK(f(1.41421) == Approx(2.24263).epsilon(1e-4));
// }
//
// TEST_CASE("Interp2D throws on extrapolation")
// {
//   constexpr const std::size_t numpoints = 2;
//   std::array<double, numpoints> const xs = {-5., 5.};
//   auto fcn = [](double x) { return -4 * x; };
//   std::array<double, numpoints> const ys = y3_cluster::transform(xs, fcn);
//   Interp2D f{xs, ys};
//   CHECK_THROWS_AS(f(-10), std::domain_error);
//   CHECK_NOTHROW(f(-5.0));
//   CHECK_NOTHROW(f(5.0));
//   CHECK_THROWS_AS(f(10), std::domain_error);
// }
