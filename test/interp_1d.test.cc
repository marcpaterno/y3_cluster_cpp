#include "catch2/catch.hpp"
#include "test/interp_1d.hh"
#include "test/transform.hh"

using y3_cluster::Interp1D;

TEST_CASE("interp_1d works", "[interpolation][1d]")
{
  constexpr const int numpoints = 9;
  std::array<double, numpoints> const xs = {1., 2., 3., 4., 5., 6, 7., 8., 9};
  auto fcn = [](double x) { return 2 * x; };
  std::array<double, numpoints> const ys = y3_cluster::transform(xs, fcn);
  Interp1D f{xs, ys};
}
