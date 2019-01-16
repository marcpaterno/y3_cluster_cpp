#include "catch2/catch.hpp"
#include "mz_power_law.hh"

#include <cmath>
#include <type_traits>

using y3_cluster::mz_power_law;

// All test answers calculated with Mathematica 11.3.0.0, using
//  f[m_, z_, B_, C_] := 1 * m^B * (1 + z)^C

TEST_CASE("has desired traits")
{
  using me = mz_power_law;
  static_assert(not std::is_default_constructible_v<me>);
  static_assert(std::is_nothrow_constructible_v<me, double, double, double>);
  static_assert(std::is_nothrow_copy_constructible_v<me>);
  static_assert(not std::is_aggregate_v<me>);
  static_assert(std::is_nothrow_invocable_r_v<double, me, double, double>);
}

TEST_CASE("works with mass = 0 and z = 0")
{
  mz_power_law f(1., 0.5, 0.5);
  CHECK(f(std::log(0.0), 1.0) == 0.0);
  CHECK(f(std::log(1.0), 0.0) == 1.0);
  CHECK(f(std::log(0.0), 0.0) == 0.0);
}

TEST_CASE("works with large mass")
{
  mz_power_law f(1, 0.1, 0.01);
  CHECK(f(std::log(1.e17), 1.0) == Approx(50.4673).epsilon(0.001));
  CHECK(f(std::log(1.e17), 10.0) == Approx(51.335).epsilon(0.001));
}

TEST_CASE("works with zero B or C")
{
  mz_power_law f{1.0, 0.0, 0.0};
  CHECK(f(std::log(1.e17), 2.0) == 1.0);
  mz_power_law g{1.0, 1.0, 0.0};
  CHECK(g(std::log(1.e13), 2.0) == Approx(1.e13).epsilon(0.001));
  mz_power_law h{1.0, 0.0, 1.0};
  CHECK(h(std::log(1.e13), 2.0) == Approx(3.0).epsilon(0.001));
}
