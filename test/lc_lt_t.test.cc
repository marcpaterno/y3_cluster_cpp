#include "catch2/catch.hpp"
#include "models/lc_lt_t.hh"
#include "utils/make_ifstream.hh"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using y3_cluster::LC_LT_t;
TEST_CASE("Lc_Lt_t works")
{
  std::ifstream infile =
    y3_cluster::make_ifstream("data/test_lc_lt_t_SDSS.txt");
  // Use REQUIRE for immediate failure if we can't open the file.
  REQUIRE(infile.good());
  std::vector<double> zs;
  std::vector<double> lts;
  std::vector<double> lcs;
  std::vector<double> ys;
  while (infile) {
    // We aren't bothering to test that the reading worked, because we're
    // careful to make sure the data file is not mal-formed when we write the
    // test.
    double z, lt, lc, y;
    infile >> z >> lt >> lc >> y;
    zs.push_back(z);
    lts.push_back(lt);
    lcs.push_back(lc);
    ys.push_back(y);
  }
  // If the file is well-formed, we have the same number of z-values as
  // y(z)-values.
  REQUIRE(zs.size() == lts.size());
  REQUIRE(lts.size() == lcs.size());
  REQUIRE(lcs.size() == ys.size());

  LC_LT_t lc_lt;

  std::ofstream out = y3_cluster::make_ofstream("validation-data/lc_lt_test.out");
  REQUIRE(out.good());
  out << std::setw(16);
  out << std::setprecision(16);
  out << "lambda_cen\tlambda_true\tprobtrue\tprobtest\n";
  for (std::size_t i = 0, sy = ys.size(); i != sy; ++i) {
    double const fz = lc_lt(lcs[i], lts[i], zs[i]);
    double constexpr epsrel = 1.0e-6;
    double constexpr epsabs = 1.0e-12;
    CHECK(fz == Approx(ys[i]).epsilon(epsrel).margin(epsabs));
    out << lcs[i] << '\t' << lts[i] << '\t' << ys[i] << '\t' << fz << '\n';
  }
}
