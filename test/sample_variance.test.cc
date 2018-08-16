#include "catch2/catch.hpp"
#include "sample_variance.hh"

using y3_cluster::SampleVariance_t,
      y3_cluster::EZ,
      y3_cluster::DV_DO_DZ_t,
      y3_cluster::IntegrationRange;

TEST_CASE("Test the sample variance calculator")
{
    auto const zz = read_vector("z_da_test.txt");
    // da_arr in h inverse Mpc
    auto const da_arr = read_vector("d_a_test.txt");
    REQUIRE(zz.size() == da_arr.size());
    auto da_f = std::make_shared<y3_cluster::Interp1D const>(zz, da_arr);

    EZ ez{0.3, 0.7, 0};
    SampleVariance_t sv(ez, {da_f, ez, 0.7}, {{{0.1, 0.3}, {0.3, 0.5}, {0.5, 0.7}}});
}
