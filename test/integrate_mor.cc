#define CATCH_CONFIG_RUNNER
#include "catch2/catch.hpp"
#include "test/probability_tests.hh"

using y3_cluster::IntegrationRange,
      y3_cluster::test_integrate_mor;
int main(int argc, char **argv)
{
    double lt_min = 1.0;
    double lt_max = 100.0;
    double M_min = 5e11;
    double M_max = 1e17;
    double zt_min = 0.1;
    double zt_max = 0.3;

    int lnM_width = 2;
    int zt_width = 2;

    int maxeval = 9999999;
    bool help = false;

    using namespace Catch::clara;
    auto args_parser = Opt(lt_min, "ltmin")
                          ["--ltmin"]
                          ("Lower bound of \\lambda_t values to try")
                     | Opt(lt_max, "ltmax")
                          ["--ltmax"]
                          ("Upper bound of \\lambda_t values to try")
                     | Opt(M_min, "Mmin")
                          ["--Mmin"]
                          ("Lower bound of mass values to try. NOTE: Will automatically convert to ln(M) to pass to mor.")
                     | Opt(M_max, "Mmax")
                          ["--Mmax"]
                          ("Upper bound of mass values to try. NOTE: Will automatically convert to ln(M) to pass to mor.")
                     | Opt(zt_min, "ztmin")
                          ["--ztmin"]
                          ("Lower bound of z_t values to try")
                     | Opt(zt_max, "ztmax")
                          ["--ztmax"]
                          ("Upper bound of z_t values to try")
                     | Opt(lnM_width, "lnM width")
                          ["--lnM-width"]
                          ("Number of distinct mass values to try")
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

    IntegrationRange lt_ir{lt_min, lt_max};
    IntegrationRange lnM_ir{std::log(M_min), std::log(M_max)};
    IntegrationRange zt_ir{zt_min, zt_max};

    cubacpp::Vegas v;
    v.maxeval = maxeval;

    test_integrate_mor(v, lt_ir, lnM_ir, zt_ir,
                       lnM_width, zt_width,
                       1.0e-4, 1.0e-12,
                       // print = true, test = false
                       true, false);
}
