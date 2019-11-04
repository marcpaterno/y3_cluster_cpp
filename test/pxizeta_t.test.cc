#include "catch2/catch.hpp"
#include "models/sptxdes/pxizeta_t.hh"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using y3_cluster::PXiZeta_t;
TEST_CASE("PXiZeta_t works")
{
  std::ifstream infile{"../data/test_pxizeta_t.txt"};

  // Test that we can open the file
  REQUIRE(infile.good());

  // Initialize vectors to hold truth
  std::vector<double> xi, zeta, prob_true;
  std::string headerline;

  // Read in the truth
  std::getline(infile, headerline);
  std::getline(infile, headerline);
  while (infile) {
    double onexi, onezeta, oneprob;
    infile >> onexi >> onezeta >> oneprob;
    xi.push_back(onexi);
    zeta.push_back(onezeta);
    prob_true.push_back(oneprob);
  }

  // Ensure that we have the some number of values in each column
  REQUIRE(xi.size() == zeta.size());
  REQUIRE(zeta.size() == prob_true.size());

  // Create an instance of the PXiZeta_t class
  PXiZeta_t pxizeta;

  // Setup output stream for validation
  std::ofstream out{"../data/pxizeta_test.out"};
  out << std::setw(16);
  out << std::setprecision(16);
  out << "xi\tzeta\tprobtest\n";

  // Test the model against truth
  double constexpr epsrel = 1.0e-6;
  double constexpr epsabs = 1.0e-12;
  for (std::size_t i = 0, len = prob_true.size(); i != len; ++i) {
    double const prob_test = pxizeta(xi[i], zeta[i]);
    CHECK(prob_test == Approx(prob_true[i]).epsilon(epsrel).margin(epsabs));
    out << xi[i] << '\t' << zeta[i] << '\t' << prob_test << '\n';
  }
}
