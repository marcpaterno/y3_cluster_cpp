#include "catch2/catch.hpp"
#include "test/omega_z_sdss.hh"

#include <fstream>

using y3_cluster::OMEGA_Z_SDSS;

TEST_CASE("omega_z_sdss works")
{
  std::ifstream infile {"omega_z_sdss.txt"};
  // Use REQUIRE for immediate failure if we can't open the file.
  REQUIRE(infile.good());
  std::vector<double> zs;
  std::vector<double> ys;
  while (infile)
  {
    // We aren't bothering to test that the reading worked, because we're
    // careful to make sure the data file is not mal-formed when we write the
    // test.
    double z, y;
    infile >> z >> y;
    zs.push_back(z);
    ys.push_back(y);
  }
  // If the file is well-formed, we have the same number of z-values as
  // y(z)-values.
  REQUIRE(zs.size() == ys.size());

  OMEGA_Z_SDSS omega;

  for (std::size_t i = 0, sz = zs.size(); i != sz; ++i)
  {
    double const fz = omega(zs[i]);
    double constexpr epsrel = 1.0e-6;
    CHECK(fz == Approx(ys[i]).epsilon(epsrel));
  }
}
