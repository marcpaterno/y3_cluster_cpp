#define CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"
#include "test/probability_tests.hh"

using y3_cluster::IntegrationRange,
      y3_cluster::test_integrate_lo_lc;
int main(int argc, char **argv)
{
    double lo_min = 1.0;
    double lo_max = 200.0;
    double lc_min = 1.0;
    double lc_max = 100.0;
    double R_min = 0.01;
    double R_max = 1.0;

    int lc_width = 2;
    int R_width = 2;

    int maxeval = 9999999;
    bool help = false;

    using namespace Catch::clara;
    auto args_parser = Opt(lo_min, "lomin")
                          ["--lomin"]
                          ("Lower bound of \\lambda_o values to try")
                     | Opt(lo_max, "lomax")
                          ["--lomax"]
                          ("Upper bound of \\lambda_o values to try")
                     | Opt(lc_min, "lcmin")
                          ["--lcmin"]
                          ("Lower bound of \\lambda_c values to try")
                     | Opt(lc_max, "lcmax")
                          ["--lcmax"]
                          ("Upper bound of \\lambda_c values to try")
                     | Opt(R_min, "R min")
                          ["--rmin"]
                          ("Lower bound of R values to try")
                     | Opt(R_max, "R max")
                          ["--rmax"]
                          ("Upper bound of R values to try")
                     | Opt(lc_width, "lc width")
                          ["--lc-width"]
                          ("Number of distinct \\lambda_c values to try")
                     | Opt(R_width, "R width")
                          ["--r-width"]
                          ("Number of distinct R values to try")
                     | Opt(maxeval, "maxeval")
                          ["--maxeval"]
                          ("Maximum number of evaluations for integration algorithm")
                     | Help(help);

    auto result = args_parser.parse(Args(argc, argv));

    if (help || !result) {
        std::cerr << args_parser;
        return 0;
    }

  IntegrationRange lo_ir{lo_min, lo_max};
  IntegrationRange lc_ir{lc_min, lc_max};
  IntegrationRange R_ir{R_min, R_max};

  cubacpp::Vegas v;
  v.maxeval = maxeval;

  test_integrate_lo_lc(v, lo_ir, lc_ir, R_ir,
                       lc_width, R_width,
                       1.0e-4, 1.0e-12,
                       // print = true, test = false
                       true, false);
}
