#include "catch2/catch.hpp"

#include <sigma_crit_inverse_t.hh>
#include <read_vector.hh>
#include <interp_1d.hh>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <stdexcept>

#include <iostream>

using namespace std;
using namespace y3_cluster;

TEST_CASE("Compare sigma_crit_inv against astropy", "[sigma_crit_inverse]")
{
  const auto z = read_vector("sigma_crit_inv/z.txt"),
             d_a = read_vector("sigma_crit_inv/d_a.txt"),
             zl = read_vector("sigma_crit_inv/zl.txt"),
             zs = read_vector("sigma_crit_inv/zs.txt"),
             astropy_sci = read_vector("sigma_crit_inv/sigma_crit_inverse.txt");
  shared_ptr<Interp1D const> da = make_shared<Interp1D const>(z, d_a);
  sigma_crit_inv sci{da};
  auto k = 0u;
  for (auto j = 0u; j < zs.size(); j++) {
    for (auto i = 0u; i < zl.size(); i++) {
      // Pop the first item from astropy_sci
      const auto expected = astropy_sci[k++];
      CHECK(sci(zl[i], zs[j]) == Approx(expected).epsilon(1e-30).margin(1e-20));
    }
  }
}

/*
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
*/
