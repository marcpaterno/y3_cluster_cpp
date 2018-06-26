#include "catch2/catch.hpp"
#include "test/interp_1d.hh"
#include "test/dv_do_dz_t.hh"
#include "test/ez.hh"
#include "test/read_vector.hh"
#include <fstream>

using y3_cluster::DV_DO_DZ_t;
using y3_cluster::Interp1D;

TEST_CASE("dv_do_dz_t works")
{
  std::ifstream infile {"test_dv_do_dz.txt"};
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

  auto identity = [](double x) { return x; };  
  auto const zz = read_vector("z_da.txt", identity);
  // da_arr in h inverse Mpc
  auto const da_arr = read_vector("d_a.txt", identity);
  REQUIRE(zz.size() == da_arr.size());
  auto da_f = std::make_shared<Interp1D const>(zz, da_arr);
  y3_cluster::DV_DO_DZ_t dvdodz(da_f, y3_cluster::EZ(0.3, 0.7, 0), 0.7);

  
  for (std::size_t i = 0, sz = zs.size(); i != sz; ++i)
  {
    double const fz = dvdodz(zs[i]);
    double constexpr epsrel = 1.0e-2;
    CHECK(fz == Approx(ys[i]).epsilon(epsrel));
  }
}
