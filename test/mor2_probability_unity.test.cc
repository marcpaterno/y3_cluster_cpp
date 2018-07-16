#include "test/probability_tests.hh"

using y3_cluster::IntegrationRange,
      y3_cluster::test_integrate_mor2;
TEST_CASE("MOR_t2 probability weighted correctly")
{
  // Lc_lt calculates $P(\lambda^{cen} | \lambda^{true}, z^{true})$ - so,
  // integral of Lc_lt over $\lambda^{cen}$ = [0, +inf] should be unity
  // Since total probability = 1
  // In practice - $\lambda^{cen}$ from 1 to 200 should be good enough
  IntegrationRange lt_ir{0.0, 1000.0};
  IntegrationRange lnM_ir{std::log(5e13), std::log(1e15)};
  IntegrationRange zt_ir{0.1, 0.3};

  cubacpp::Vegas v;
  v.maxeval = 9999999;

  test_integrate_mor2(v, lt_ir, lnM_ir, zt_ir, 2, 2);
}
