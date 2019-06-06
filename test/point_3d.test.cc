#include "catch2/catch.hpp"
#include "utils/point_3d.hh"

#include <cmath>
#include <limits>

using fpsupport::is_equivalent;
using y3_cluster::Point3D;
using y3_cluster::Point3DLess;

// avoiding too much typing...
auto infinity = std::numeric_limits<double>::infinity;

TEST_CASE("equivalence stuff")
{
  double const a = 1.e15;
  double const b = std::nextafter(a, infinity());
  CHECK(not(a == b));
  CHECK(a < b);
  double const absdiff = b - a;
  double const reldiff = absdiff / std::max(a, b);
  CHECK(is_equivalent(a, b, reldiff, absdiff));

  double const c = std::nextafter(b, infinity());
  CHECK(not is_equivalent(a, c, reldiff, absdiff));
}

TEST_CASE("sorting stuff")
{
  Point3DLess p3dless{1.e-6, 1.e-20, 1.e-6, 1.e-20};

  // Construct so that a < b < c.
  Point3D const a{2.25, -5.0, -100.};
  Point3D const b{1.5, 2.0, -200.};
  Point3D const c{1.25, 4.0, -300.};

  SECTION("ordering")
  {
    CHECK(p3dless(a, b));
    CHECK(p3dless(b, c));
    CHECK(p3dless(a, c));
    CHECK(not p3dless(a, a));
    CHECK(not p3dless(b, b));
    CHECK(not p3dless(c, c));
    CHECK(not p3dless(b, a));
    CHECK(not p3dless(c, b));
    CHECK(not p3dless(c, a));
  }

  SECTION("sorting")
  {
    std::vector<Point3D> pts{c, b, a};
    CHECK(not std::is_sorted(pts.begin(), pts.end(), p3dless));
    std::sort(pts.begin(), pts.end(), p3dless);
    CHECK(std::is_sorted(pts.begin(), pts.end(), p3dless));
  }
}
