#define CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"
#include "test/probability_tests.hh"

using y3_cluster::IntegrationRange,
      y3_cluster::test_integrate_lo_lc;
int main(void)
{
  IntegrationRange lo_ir{1.0, 200.0};
  IntegrationRange lc_ir{1.0, 100.0};
  IntegrationRange R_ir{0.01, 1.0};

  cubacpp::Vegas v;
  v.maxeval = 9999999;

  test_integrate_lo_lc(v, lo_ir, lc_ir, R_ir, 2, 2,
                       1.0e-4, 1.0e-12,
                       // print = true, test = false
                       true, false);
}
