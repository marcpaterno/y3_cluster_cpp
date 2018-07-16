#include "catch2/catch.hpp"
#include "hmf_t.hh"
#include "interp_2d.hh"
#include "read_vector.hh"
#include <fstream>

using y3_cluster::HMF_t;
using y3_cluster::Interp2D;

TEST_CASE("mass function works")
{
  std::ifstream infile {"test_hmf_z0_z03.dat"};
  // Use REQUIRE for immediate failure if we can't open the file.
  REQUIRE(infile.good());
  std::vector<double> zs;
  std::vector<double> ms;
  std::vector<double> ys;
  while (infile)
  {
    // We aren't bothering to test that the reading worked, because we're
    // careful to make sure the data file is not mal-formed when we write the
    // test.
    double z, m, y;
    infile >> z >> m >> y;
    zs.push_back(z);
    ms.push_back( m*std::log(10) );
    ys.push_back(y*std::pow(10, m));
  }
  // If the file is well-formed, we have the same number of z-values as
  // y(z)-values.
  REQUIRE(zs.size() == ys.size());
  REQUIRE(zs.size() == ms.size());

  const double Omega_M = 1.87518978e-01;
  auto identity = [](double x) { return x; };
  auto log_omega_m = [Omega_M](double x) { return std::log(x*Omega_M); };
  auto mh = read_vector("m_h_matteo.txt", log_omega_m);
  auto const zz = read_vector("z_matteo.txt", identity);
  auto const dndlnmh = read_vector("dndlnmh_matteo.txt", identity);
  auto p1 = std::make_shared<Interp2D const>(mh, zz, dndlnmh);
  //HMF_t hmf(p1, 4.50732047e-02, 1.01958078e+00);
  HMF_t hmf(p1, 0.0, 1.0);

  for (std::size_t i = 0, sz = zs.size(); i != sz; ++i)
  {
    double const fz = hmf(ms[i], zs[i]);
    double constexpr epsrel = 1.0e-3;
    CHECK(fz == Approx(ys[i]).epsilon(epsrel));
  }
}
