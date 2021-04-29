#include "catch2/catch.hpp"
#include "utils/interp_2d.hh"
#include "utils/transform.hh"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <stdexcept>

#include <iostream>

#include "gsl/gsl_errno.h"

using y3_cluster::Interp2D;
using y3_cluster::Point3D;
auto infinity = std::numeric_limits<double>::infinity;

template <std::size_t M, std::size_t N>
std::size_t
nentries(Interp2D::matrix<M, N> const&)
{
  return M * N;
}

TEST_CASE("clamp interface works")
{
  constexpr std::size_t nx = 3;
  constexpr std::size_t ny = 2;
  std::array<double, nx> const xs = {1., 2., 3.};
  std::array<double, ny> const ys = {4., 5.};
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

  SECTION("interpolation works") { CHECK(f(2.5, 4.5) == 56.75); }
  SECTION("extrapolation gets clamped")
  {
    CHECK(f.clamp(0, 4.5) == f(1, 4.5)); // to the left
    CHECK(f.clamp(4, 4.5) == f(3, 4.5)); // to the right
    CHECK(f.clamp(2, 3) == f(2, 4));     // below
    CHECK(f.clamp(2, 5.5) == f(2, 5));   // above
    CHECK(f.clamp(0, 0) == f(1., 4.));   // below left
    CHECK(f.clamp(4, 3) == f(3, 4));     // below right
    CHECK(f.clamp(0, 6) == f(1, 5));     // above left
    CHECK(f.clamp(4, 6) == f(3, 5));     // above right
  }
}

TEST_CASE("Interp2D exact at knots", "[interpolation][2d]")
{
  constexpr std::size_t nx = 3;
  constexpr std::size_t ny = 2;
  std::array<double, nx> const xs = {1., 2., 3.};
  std::array<double, ny> const ys = {4., 5.};
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
  gsl_set_error_handler_off();
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
  /*
   * FIXME: Add these tests back when we replace the Interp2D bounds checking.
  CHECK_THROWS_AS(f(0.0, 3.0) == 4.0, std::domain_error);
  CHECK_THROWS_AS(f(5.0, 3.0) == 14.0, std::domain_error);
  // Extrapolate, only in y.
  CHECK_THROWS_AS(f(2.5, -1.0) == -3.0, std::domain_error);
  CHECK_THROWS_AS(f(2.5, 5.0) == 15.0, std::domain_error);
  // Extrapolate in both x and y.
  CHECK_THROWS_AS(f(0.0, 0.0) == -5.0, std::domain_error);
  CHECK_THROWS_AS(f(5.0, 5.0) == 20.0, std::domain_error);
  */
}

TEST_CASE("construction from Point3Ds")
{
  SECTION("from the minimal grid")
  {
    std::array<double, 2> const xs{1, 2};
    std::array<double, 2> const ys{4, 5};
    std::vector<Point3D> points;
    auto fxy = [](double x, double y) { return 2 * x + 3 * y - 5; };

    for (auto x : xs) {
      for (auto y : ys) {
        Point3D p{x, y, fxy(x, y)};
        points.push_back(p);
      }
    }
    CHECK(points.size() == 4);

    Interp2D::matrix<2, 2> zs;
    for (std::size_t i = 0; i != 2; ++i) {
      double x = xs[i];
      for (std::size_t j = 0; j != 2; ++j) {
        double y = ys[j];
        zs[i][j] = fxy(x, y);
        CHECK(zs[i][j] == fxy(x, y));
      }
    }
    Interp2D f2{xs, ys, zs};

    CHECK_NOTHROW(Interp2D(points));
    Interp2D f(points);
    for (std::size_t i = 0; i != xs.size(); ++i) {
      double const x = xs[i];
      for (std::size_t j = 0; j != ys.size(); ++j) {
        double const y = ys[j];
        CHECK(f2(x, y) == fxy(x, y));
        CHECK(f(x, y) == fxy(x, y));
      }
    }
  }
  SECTION("from a correct grid")
  {
    std::array<double, 3> const xs{2.5, 3.25, 5.8};
    std::array<double, 2> const ys{-1.3, 2.4};
    auto fcn = [](double x, double y) {
      return std::sin(x) + 2 * x * y + y * y;
    };
    std::vector<Point3D> points;
    points.reserve(xs.size() * ys.size());
    for (std::size_t i = 0; i != xs.size(); ++i) {
      double const x = xs[i];
      for (std::size_t j = 0; j != ys.size(); ++j) {
        double const y = ys[j];
        Point3D const newguy{x, y, fcn(x, y)};
        points.push_back(newguy);
      }
    }
    CHECK(points.size() == xs.size() * ys.size());
    CHECK_NOTHROW(Interp2D(points));
    std::reverse(points.begin(), points.end());
    std::mt19937 g(271); // not worried about a very good seed.
    std::shuffle(points.begin(), points.end(), g);
    CHECK_NOTHROW(Interp2D(points));
  }
}

TEST_CASE("detect defective grids")
{
  SECTION("no points provided")
  {
    // First lvalue reference...
    std::vector<Point3D> const no_points;
    CHECK_THROWS_AS(Interp2D(no_points), std::domain_error);
    // Now rvalue reference...
    CHECK_THROWS_AS(Interp2D(std::vector<Point3D>()), std::domain_error);
  }
  SECTION("missing an entry")
  {
    std::array<double, 3> const xs{2.5, 3.25, 5.8};
    std::array<double, 2> const ys{-1.3, 2.4};
    auto fcn = [](double x, double y) {
      return std::sin(x) + 2 * x * y + y * y;
    };
    std::vector<Point3D> points;
    points.reserve(xs.size() * ys.size());
    for (std::size_t i = 0; i != xs.size(); ++i) {
      double const x = xs[i];
      for (std::size_t j = 0; j != ys.size(); ++j) {
        if (i == 1 && j == 1) continue; // skip one point
        double const y = ys[j];
        Point3D const newguy{x, y, fcn(x, y)};
        points.push_back(newguy);
      }
    }
    CHECK_THROWS_AS(Interp2D(points), std::domain_error);
  }
  SECTION("slightly wobbly axis values are ok")
  {
    std::array<double, 3> xs{2.5, 3.25, 5.8};
    std::array<double, 2> ys{-1.3, 2.4};
    auto fcn = [](double x, double y) {
      return std::sin(x) + 2 * x * y + y * y;
    };
    std::vector<Point3D> points;
    points.reserve(xs.size() * ys.size());
    for (std::size_t i = 0; i != xs.size(); ++i) {
      double const x = xs[i];
      for (std::size_t j = 0; j != ys.size(); ++j) {
        double const y = ys[j];
        Point3D const newguy{x, y, fcn(x, y)};
        points.push_back(newguy);
      }
    }
    CHECK(points.size() > 3); // sanity test before tweaking points.
    double const xoriginal = points[3][0];
    points[3][0] = std::nextafter(xoriginal, infinity());
    double const yoriginal = points[2][1];
    points[2][1] = std::nextafter(yoriginal, -infinity());
    CHECK_NOTHROW(Interp2D(points));
  }
  SECTION("badly wobbly axis values are no good")
  {
    std::array<double, 3> xs{2.5, 3.25, 5.8};
    std::array<double, 2> ys{-1.3, 2.4};
    auto fcn = [](double x, double y) {
      return std::sin(x) + 2 * x * y + y * y;
    };
    std::vector<Point3D> points;
    points.reserve(xs.size() * ys.size());
    for (std::size_t i = 0; i != xs.size(); ++i) {
      double const x = xs[i];
      for (std::size_t j = 0; j != ys.size(); ++j) {
        double const y = ys[j];
        Point3D const newguy{x, y, fcn(x, y)};
        points.push_back(newguy);
      }
    }
    CHECK(points.size() > 3); // sanity test before tweaking points.
    points[3][0] *= (1. - 2.e-6);
    points[2][1] *= (1. + 2.e-6);
    CHECK_THROWS_AS(Interp2D(points), std::domain_error);
  }
}
