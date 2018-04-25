#include "catch2/catch.hpp"
#include "test/fpsupport.hh"
#include "test/point_3d.hh"

#include <cmath>
#include <limits>

#include <iomanip>
#include <iostream>

using fpsupport::icky;

// avoiding too much typing...
auto infinity = std::numeric_limits<double>::infinity;
auto NaN = std::numeric_limits<double>::quiet_NaN;
auto tiny = std::numeric_limits<double>::denorm_min;
auto small = std::numeric_limits<double>::min;

TEST_CASE("icky doubles")
{
  auto const bad = 0.0 / 0.0;
  CHECK(std::isnan(bad));
  CHECK(std::numeric_limits<double>::has_quiet_NaN);
  auto const a = std::fpclassify(NaN());
  CHECK(std::isnan(NaN()));
  CHECK(a == FP_NAN);
  CHECK(not icky(0.0));
  CHECK(not icky(-1.25));
  CHECK(not icky(1.23e204));
  CHECK(not icky(tiny()));
  CHECK(not icky(small()));
  CHECK(icky(NaN()));
  CHECK(icky(-infinity()));
  CHECK(icky(infinity()));
}

TEST_CASE("subnormals behave as expected")
{
  CHECK(not std::isnormal(tiny()));
  CHECK(not std::isnormal(0.0));
  CHECK(std::isnormal(small()));
  CHECK(std::fpclassify(0.0) == FP_ZERO);
  CHECK(std::fpclassify(tiny()) != FP_ZERO);
  CHECK(std::fpclassify(small()) != FP_ZERO);
  double const a = 2.0 * tiny();
  CHECK(a != tiny());
  CHECK(a > tiny());
  CHECK(tiny() < a);
  double const b = tiny() + tiny();
  CHECK(b != tiny());
  CHECK(b > tiny());
  CHECK(tiny() < b);
}

TEST_CASE("squashing works")
{
  y3_cluster::Point3D p1{tiny(), tiny(), tiny()};
  for (auto x : p1) {
    CHECK(not std::isnormal(x));
    CHECK(std::fpclassify(x) != FP_ZERO);
  }
  y3_cluster::squash_subnormals(p1);
  CHECK(std::fpclassify(p1[0]) == FP_ZERO);
  CHECK(std::fpclassify(p1[1]) == FP_ZERO);
  CHECK(std::fpclassify(p1[2]) != FP_ZERO);

  y3_cluster::Point3D p2{-tiny(), -tiny(), -tiny()};
  for (auto x : p2) {
    CHECK(not std::isnormal(x));
    CHECK(std::fpclassify(x) != FP_ZERO);
  }
  y3_cluster::squash_subnormals(p2);
  CHECK(std::fpclassify(p2[0]) == FP_ZERO);
  CHECK(std::fpclassify(p2[1]) == FP_ZERO);
  CHECK(std::fpclassify(p2[2]) != FP_ZERO);
}
