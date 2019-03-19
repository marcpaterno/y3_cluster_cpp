#include "catch2/catch.hpp"
#include "primitives.hh"

using y3_cluster::integer_pow,
      y3_cluster::pi;

TEST_CASE("integer_pow works")
{
  std::vector<double> values{{0.1, 1.0, pi(), 127.324, -195.3243}};
  for (const auto v : values) {
      double running_up = 1.0, running_down = 1.0;
      for (auto i = 0; i < 20; i++) {
          CHECK(integer_pow(v, i) == Approx(running_up).epsilon(1e-10));
          CHECK(integer_pow(v, -i) == Approx(running_down).epsilon(1e-10));
          running_up *= v;
          running_down /= v;
      }
  }
}
