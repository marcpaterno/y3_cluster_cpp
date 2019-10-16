#include "catch2/catch.hpp"
#include "utils/polynomial.hh"

using y3_cluster::polynomial;

TEST_CASE("Polynomial implementation")
{
  polynomial<3> cubic{{2, -1, 3}};
  for (auto i = 0u; i < 100; i++) {
    double x = (i - 50) / (50.0 / 3.0);
    CHECK(cubic(x) == Approx(2 * x * x - x + 3).epsilon(1e-8));
  }

  polynomial<4> quartic{{-5, -2, 7, 2.3}};
  for (auto i = 0u; i < 100; i++) {
    double x = (i - 50) / (50.0 / 3.0);
    CHECK(quartic(x) == Approx(((-5 * x - 2) * x + 7) * x + 2.3).epsilon(1e-8));
  }

  polynomial<5> quintic{{4, 2.75, -1, 0, -11.3}};
  for (auto i = 0u; i < 100; i++) {
    double x = (i - 50) / (50.0 / 3.0);
    CHECK(quintic(x) ==
          Approx(((4 * x + 2.75) * x - 1) * x * x - 11.3).epsilon(1e-8));
  }
}
