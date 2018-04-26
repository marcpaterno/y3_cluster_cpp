#include "catch2/catch.hpp"
#include "test/interp_2d.hh"
#include "test/transform.hh"

#include <cmath>
#include <cstddef>
#include <stdexcept>

using y3_cluster::Interp2D;

template <std::size_t M, std::size_t N>
std::size_t
nentries(Interp2D::matrix<M, N> const&)
{
  return M * N;
}

TEST_CASE("Interp2D exact at knots", "[interpolation][2d]")
{
  constexpr std::size_t nx = 3;
  constexpr std::size_t ny = 2;
  std::array<double, nx> const xs = {1., 2., 3.};
  std::array<double, ny> const ys = {1., 2.};
  auto fxy = [](double x, double y) { return 3 * x * y + 2 * x + 4 * y; };
  Interp2D::matrix<nx, ny> zs;
  CHECK(nentries(zs) == nx * ny);
  CHECK(sizeof(zs) == nx * ny * sizeof(double));
  for (std::size_t i = 0; i != nx; ++i) {
    double x = xs[i];
    for (std::size_t j = 0; j != ny; ++j) {
      double y = ys[j];
      zs[i][j] = fxy(x, y);
      CHECK(zs[i][j] == fxy(x, y));
    }
  }

  Interp2D f{xs, ys, zs};
  for (std::size_t i = 0; i != nx; ++i) {
    double x = xs[i];
    for (std::size_t j = 0; j != ny; ++j) {
      double y = ys[j];
      CHECK(zs[i][j] == fxy(x, y));
      CHECK(zs[i][j] == f(x, y));
    }
  }
}

TEST_CASE("Interp2D on bilinear")
{
  constexpr std::size_t nx = 3;
  constexpr std::size_t ny = 4;
  std::array<double, nx> const xs = {1., 2., 3.};
  std::array<double, ny> const ys = {1., 2., 3., 4.};

  auto fxy = [](double x, double y) { return 2 * x + 3 * y - 5; };
  Interp2D::matrix<nx, ny> zs;

  for (std::size_t i = 0; i != nx; ++i) {
    double x = xs[i];
    for (std::size_t j = 0; j != ny; ++j) {
      double y = ys[j];
      zs[i][j] = fxy(x, y);
      CHECK(zs[i][j] == fxy(x, y));
    }
  }

  Interp2D f{xs, ys, zs};
  // Pure interpolation.
  CHECK(f(2.5, 1.5) == 4.5);
  // Extrapolate, only in x.
  CHECK(f(0.0, 3.0) == 4.0);
  CHECK(f(5.0, 3.0) == 14.0);
  // Extrapolate, only in y.
  CHECK(f(2.5, -1.0) == -3.0);
  CHECK(f(2.5, 5.0) == 15.0);
  // Extrapolate in both x and y.
  CHECK(f(0.0, 0.0) == -5.0);
  CHECK(f(5.0, 5.0) == 20.0);
}
