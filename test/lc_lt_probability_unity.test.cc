#include "catch2/catch.hpp"
#include "integration_range.hh"
#include "cubacpp/cubacpp.hh"
#include "test/lc_lt_t.hh"
#include "test/probability_tests.hh"

#include <fstream>
#include <iostream>
#include <string>

using y3_cluster::LC_LT_t,
      y3_cluster::IntegrationRange,
      y3_cluster::test_integrate_lc_lt;
TEST_CASE("Lc_Lt_t probability weighted correctly")
{
  // Lc_lt calculates $P(\lambda^{cen} | \lambda^{true}, z^{true})$ - so,
  // integral of Lc_lt over $\lambda^{cen}$ = [0, +inf] should be unity
  // Since total probability = 1
  // In practice - $\lambda^{cen}$ from 1 to 200 should be good enough
  IntegrationRange lc_ir{1.0, 200.0};
  IntegrationRange lt_ir{1.0, 100.0};
  IntegrationRange zt_ir{0.1, 0.3};

  cubacpp::Vegas v;
  v.maxeval = 9999999;

  test_integrate_lc_lt(v, lc_ir, lt_ir, zt_ir, 2, 2);
}
