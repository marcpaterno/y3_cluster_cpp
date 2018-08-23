#include "catch2/catch.hpp"
#include "sample_variance.hh"

using y3_cluster::SampleVariance_t,
      y3_cluster::IntegrationRange;

TEST_CASE("Test the sample variance calculator")
{
    cubacores(0, 0);
    const double h = 0.771358152;
    auto const zz = read_vector("z_da_test.txt");
    // da_arr in h inverse Mpc
    auto const da_arr = read_vector("d_a_test.txt");
    REQUIRE(da_arr.size() == zz.size());
    auto da_f = std::make_shared<y3_cluster::Interp1D const>(zz, da_arr);

    SampleVariance_t sv({{{0.1, 0.3}}}, h);

    auto matrix = sv.compute();
    REQUIRE(matrix.size() == 1);
    for (const auto& row : matrix)
        REQUIRE(row.size() == 1);
}
