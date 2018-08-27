#include "test/probability_tests.hh"


using y3_cluster::IntegrationRange,
      y3_cluster::test_integrate_roffset;
TEST_CASE("roffset probability weighted correctly")
{
    IntegrationRange R_ir{0.0, 1.0};

    cubacpp::CQUAD v;

    test_integrate_roffset(v, R_ir, 0.2);
}
