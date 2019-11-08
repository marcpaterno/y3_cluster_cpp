#include "catch2/catch.hpp"
#include "models/sptxdes/int_pxizeta.hh"

#include <fstream>
#include <iostream>
#include <string>

using y3_cluster::INT_PXIZETA;

TEST_CASE("Testing Integral over xi of P(xi|zeta)")
{
  std::ifstream infile{"../data/test_int_pxizeta.txt"};

  // Test that we can open the file
  REQUIRE(infile.good());

  // Initialize vectors to hold the truth
  std::vector<double> zeta, truth_integral;
  std::string headerline;

  // Read in the truth
  std::getline(infile, headerline);
  std::getline(infile, headerline);
  while(infile) {
    double onezeta, onetruth;
    infile >> onezeta >> onetruth;
    zeta.push_back(onezeta);
    truth_integral.push_back(onetruth);
  }

  // Ensure we have the same number of truth values per column
  REQUIRE(zeta.size() == truth_integral.size());

  // Create an instance of the class to test
  INT_PXIZETA int_pxizeta;

  // Test against the truth
  double constexpr epsrel = 1.0e-6;
  double constexpr epsabs = 1.0e-12;
  for (std::size_t i = 0, len = zeta.size(); i != len; ++i) {
    double const test_integral = int_pxizeta(zeta[i]);
    CHECK(test_integral == Approx(truth_integral[i]).epsilon(epsrel).margin(epsabs));
  }
}
