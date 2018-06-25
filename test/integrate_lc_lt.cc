#define CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"
#include "test/probability_tests.hh"

using y3_cluster::IntegrationRange,
      y3_cluster::test_integrate_lc_lt;
int main(int argc, char **argv)
{
    double lc_min = 1.0;
    double lc_max = 200.0;
    double lt_min = 1.0;
    double lt_max = 100.0;
    double zt_min = 0.1;
    double zt_max = 0.3;

    int lt_width = 2;
    int zt_width = 2;

    int maxeval = 9999999;
    bool help = false;

    using namespace Catch::clara;
    auto args_parser = Opt(lc_min, "lcmin")
                          ["--lcmin"]
                          ("Lower bound of \\lambda_c values to try")
                     | Opt(lc_max, "lcmax")
                          ["--lcmax"]
                          ("Upper bound of \\lambda_c values to try")
                     | Opt(lt_min, "ltmin")
                          ["--ltmin"]
                          ("Lower bound of \\lambda_t values to try")
                     | Opt(lt_max, "ltmax")
                          ["--ltmax"]
                          ("Upper bound of \\lambda_t values to try")
                     | Opt(lt_min, "ztmin")
                          ["--ztmin"]
                          ("Lower bound of z_t values to try")
                     | Opt(lt_max, "ztmax")
                          ["--ztmax"]
                          ("Upper bound of z_t values to try")
                     | Opt(lt_width, "lt width")
                          ["--lt-width"]
                          ("Number of distinct \\lambda_t values to try")
                     | Opt(zt_width, "zt width")
                          ["--zt-width"]
                          ("Number of distinct z_t values to try")
                     | Opt(maxeval, "maxeval")
                          ["--maxeval"]
                          ("Maximum number of evaluations for integration algorithm")
                     | Help(help);

    auto result = args_parser.parse(Args(argc, argv));

    if (help || !result) {
        std::cerr << args_parser;
        return 0;
    }

    IntegrationRange lc_ir{lc_min, lc_max};
    IntegrationRange lt_ir{lt_min, lt_max};
    IntegrationRange zt_ir{zt_min, zt_max};

    cubacpp::Vegas v;
    v.maxeval = maxeval;

    test_integrate_lc_lt(v, lc_ir, lt_ir, zt_ir,
                         lt_width, zt_width,
                         1.0e-4, 1.0e-12,
                         // print = true, test = false
                         true, false);
}
