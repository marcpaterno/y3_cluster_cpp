#include "test/probability_tests.hh"

using y3_cluster::LO_LC_t,
      y3_cluster::IntegrationRange,
      y3_cluster::test_integrate_lo_lc;
TEST_CASE("LO_LC_t probability weighted correctly")
{
  // Lc_lt calculates $P(\lambda^{cen} | \lambda^{true}, z^{true})$ - so,
  // integral of Lc_lt over $\lambda^{cen}$ = [0, +inf] should be unity
  // Since total probability = 1
  // In practice - $\lambda^{cen}$ from 1 to 200 should be good enough
  IntegrationRange lo_ir{1.0, 200.0};
  IntegrationRange lc_ir{1.0, 100.0};
  IntegrationRange R_ir{0.01, 1.0};

  cubacpp::Vegas v;
  v.maxeval = 9999999;

  test_integrate_lo_lc(v, lo_ir, lc_ir, R_ir, 2, 2);
}
