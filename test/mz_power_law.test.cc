#include "catch2/catch.hpp"
#include "test/mz_power_law.hh"

#include <type_traits>

using y3_cluster::mz_power_law;

TEST_CASE("has desired traits")
{
  static_assert(not std::is_default_constructible_v<mz_power_law>);
  static_assert(
    std::is_nothrow_constructible_v<mz_power_law, double, double, double>);
  static_assert(std::is_nothrow_copy_constructible_v<mz_power_law>);
  static_assert(not std::is_aggregate_v<mz_power_law>);
  static_assert(
    std::is_nothrow_invocable_r_v<double, mz_power_law, double, double>);
}
