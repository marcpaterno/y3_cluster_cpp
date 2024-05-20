#include "catch2/catch.hpp"
#include "models/lo_lc_t.hh"
#include "utils/make_ifstream.hh"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using y3_cluster::LO_LC_t;
TEST_CASE("Lo_Lc_t works")
{
  std::ifstream infile = y3_cluster::make_ifstream("data/test_lolc_t.txt");
  // Use REQUIRE for immediate failure if we can't open the file.
  REQUIRE(infile.good());
  std::vector<double> los;
  std::vector<double> lcs;
  std::vector<double> Rs;
  std::vector<double> ys;
  std::string dummyline;
  std::getline(infile, dummyline);
  std::getline(infile, dummyline);
  while (infile) {
    // We aren't bothering to test that the reading worked, because we're
    // careful to make sure the data file is not mal-formed when we write the
    // test.
    double lo, lc, R, y;
    infile >> lo >> lc >> R >> y;
    los.push_back(lo);
    lcs.push_back(lc);
    Rs.push_back(R);
    ys.push_back(y);
  }
  // If the file is well-formed, we have the same number of z-values as
  // y(z)-values.
  REQUIRE(los.size() == lcs.size());
  REQUIRE(lcs.size() == Rs.size());
  REQUIRE(Rs.size() == ys.size());

  std::ofstream out = y3_cluster::make_ofstream("data/test_lolc.out");
  REQUIRE(out.good());
  out << std::setw(16);
  out << std::setprecision(16);
  out << "lambda_o\tlambda_c\tR_mis\tytrue\tytest\n";

  // No longer relevant - redefined lo_lc
  LO_LC_t lolc(1.66, 0.26, 1.43, 1.0);

  std::ofstream out2 = y3_cluster::make_ofstream("validation-data/lo_lc_test.out");
  REQUIRE(out2.good());
  out2 << std::setw(16);
  out2 << std::setprecision(16);
  out2 << "deltal\tRmis\tprobtrue\tprobtest\n";
  for (std::size_t i = 0, sy = ys.size(); i != sy; ++i) {
    double const fz = lolc(los[i], lcs[i], Rs[i]);
    double constexpr epsrel = 1.0e-6;
    double constexpr epsabs = 1.0e-12;
    CHECK((fz * lcs[i]) == Approx(ys[i]).epsilon(epsrel).margin(epsabs));
    out2 << los[i] / lcs[i] << '\t' << Rs[i] << '\t' << ys[i] << '\t'
         << fz * lcs[i] << '\n';
  }
}
