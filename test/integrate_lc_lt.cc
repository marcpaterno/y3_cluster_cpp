#define CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"
#include "test/probability_tests.hh"

using y3_cluster::IntegrationRange,
      y3_cluster::test_integrate_lc_lt;
int main(void)
{
  IntegrationRange lc_ir{1.0, 200.0};
  IntegrationRange lt_ir{1.0, 100.0};
  IntegrationRange zt_ir{0.1, 0.3};

  cubacpp::Vegas v;
  v.maxeval = 9999999;

  test_integrate_lc_lt(v, lc_ir, lt_ir, zt_ir, 2, 2,
                       1.0e-4, 1.0e-12,
                       // print = true, test = false
                       true, false);
}
