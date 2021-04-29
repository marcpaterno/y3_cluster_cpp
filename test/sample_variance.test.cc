#include <cubacpp/cubacpp.hh>

#include "catch2/catch.hpp"
#include "models/omega_z_sdss.hh"
#include "models/sample_variance.hh"

using y3_cluster::IntegrationRange, y3_cluster::OMEGA_Z_SDSS,
  y3_cluster::SampleVariance_t;

TEST_CASE("Test the sample variance calculator")
{
  cubacpp::turn_off_cuba_forking();
  const double h = 0.771358152;
  auto const zz = read_vector("z_da_test.txt");
  // da_arr in h inverse Mpc
  auto const da_arr = read_vector("d_a_test.txt");
  REQUIRE(da_arr.size() == zz.size());
  auto da_f = std::make_shared<y3_cluster::Interp1D const>(zz, da_arr);

  OMEGA_Z_SDSS omega_z;
  SampleVariance_t sv(omega_z, {{{0.1, 0.3}, {0.3, 0.5}, {0.5, 0.7}}}, h);

  auto matrix = sv.compute();
  REQUIRE(matrix.size() == 3);
  for (const auto& row : matrix) REQUIRE(row.size() == 3);
}
