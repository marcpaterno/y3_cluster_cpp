#include "catch2/catch.hpp"
#include <fstream>
#include <iomanip>
#include <nfw.hh>


TEST_CASE("Xi_NFW works")
{
  std::ifstream infile {"xi_14.5_nfw.txt"};
  // Use REQUIRE for immediate failure if we can't open the file.
  REQUIRE(infile.good());
  std::vector<double> Rs;
  std::vector<double> ys;
  while (infile)
  {
    // We aren't bothering to test that the reading worked, because we're
    // careful to make sure the data file is not mal-formed when we write the
    // test.
    double R, y;
    infile >> R >> y;
    Rs.push_back(R);
    ys.push_back(y);
    std::cout << R << y;
  }
  // If the file is well-formed, we have the same number of z-values as
  // y(z)-values.
  REQUIRE(Rs.size() == ys.size());

  nfw_dsigma nfw(5.0, 0.3);

  std::ofstream out { "../data/xi_14.5_nfw.out" };
  out << std::setw(16);
  out << std::setprecision(16);
  out << "z\tytrue\tytest\n";

  for (std::size_t i = 0, sR = Rs.size(); i != sR; ++i)
  {
    double const fz = nfw(Rs[i], pow(10.0, 14.5) );
    double constexpr epsrel = 1.0e-3;
    CHECK(fz == Approx(ys[i]).epsilon(epsrel));
    out << Rs[i] << '\t' << ys[i] << '\t' << fz << '\n';
  }
}
