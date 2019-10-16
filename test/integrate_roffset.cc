#define CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"
#include "test/probability_tests.hh"

using y3_cluster::IntegrationRange, y3_cluster::test_integrate_roffset;
int
main(int argc, char** argv)
{
  double R_min = 0.0;
  double R_max = 1.0;
  double tau = 0.2;
  int maxeval = 9999999;
  bool help = false;

  using namespace Catch::clara;
  auto args_parser =
    Opt(R_min, "R min")["--r-min"]("Lower integration bound of R") |
    Opt(R_max, "R max")["--r-max"]("Upper integration bound of R") |
    Opt(tau,
        "tau")["--tau"]("Width parameter of roffset. Only one value used") |
    Opt(maxeval, "maxeval")["--maxeval"](
      "Maximum number of evaluations for integration algorithm") |
    Help(help);

  auto result = args_parser.parse(Args(argc, argv));

  if (help || !result) {
    std::cerr << args_parser;
    return 0;
  }

  IntegrationRange R_ir{R_min, R_max};

  cubacpp::Vegas v;
  v.maxeval = maxeval;

  test_integrate_roffset(v,
                         R_ir,
                         tau,
                         1.0e-4,
                         1.0e-12,
                         // print = true, test = false
                         true,
                         false);
}
