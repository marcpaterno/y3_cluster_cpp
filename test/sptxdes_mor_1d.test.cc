#include "catch2/catch.hpp"
#include "models/sptxdes/mor_1d.hh"
#include "utils/make_ifstream.hh"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using y3_cluster::MOR_1D;
TEST_CASE("SPTxDES MOR_1D works")
{
  std::ifstream infile = y3_cluster::make_ifstream("data/test_sptxdes_mor_1d.txt");

  // Test that we can open the file
  REQUIRE(infile.good());

  // Initialize vectors to hold truth
  std::vector<double> mmin, m1, alpha, epsilon, sigintr;
  std::vector<double> zplam, lamtrue, ztrue, lnm200m, probtrue;
  std::string headerline;

  // Read in the truth
  std::getline(infile, headerline);
  std::getline(infile, headerline);
  while (infile) {
    double onemmin, onem1, onealpha, oneepsilon, onesigintr;
    double onezp, onelt, onezt, onelnm, oneprob;
    infile >> onemmin >> onem1 >> onealpha >> oneepsilon >> onesigintr;
    infile >> onezp >> onelt >> onezt >> onelnm >> oneprob;
    mmin.push_back(onemmin);
    m1.push_back(onem1);
    alpha.push_back(onealpha);
    epsilon.push_back(oneepsilon);
    sigintr.push_back(onesigintr);
    zplam.push_back(onezp);
    lamtrue.push_back(onelt);
    ztrue.push_back(onezt);
    lnm200m.push_back(onelnm);
    probtrue.push_back(oneprob);
  }

  // Ensure that we have the some number of values in each column
  REQUIRE(mmin.size() == m1.size());
  REQUIRE(mmin.size() == alpha.size());
  REQUIRE(mmin.size() == epsilon.size());
  REQUIRE(mmin.size() == sigintr.size());
  REQUIRE(mmin.size() == zplam.size());
  REQUIRE(mmin.size() == lamtrue.size());
  REQUIRE(mmin.size() == ztrue.size());
  REQUIRE(mmin.size() == lnm200m.size());
  REQUIRE(mmin.size() == probtrue.size());

  // Test the model against truth
  double constexpr epsrel = 1.0e-6;
  double constexpr epsabs = 1.0e-12;
  for (std::size_t i = 0, len = probtrue.size(); i != len; ++i) {
    MOR_1D mor(mmin[i], m1[i], alpha[i], epsilon[i], sigintr[i], zplam[i]);
    double const probtest = mor(lamtrue[i], ztrue[i], lnm200m[i]);
    CHECK(probtest == Approx(probtrue[i]).epsilon(epsrel).margin(epsabs));
  }
}
