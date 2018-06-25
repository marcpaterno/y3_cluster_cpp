#define CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"
#include "test/probability_tests.hh"

using y3_cluster::IntegrationRange,
      y3_cluster::test_integrate_roffset;
int main(void)
{
    IntegrationRange R_ir{0.0, 1.0};

    cubacpp::Vegas v;
    v.maxeval = 9999999;

    test_integrate_roffset(v, R_ir, 0.2,
                           1.0e-4, 1.0e-12,
                           // print = true, test = false
                           true, false);
}
