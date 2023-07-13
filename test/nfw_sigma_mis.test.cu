#include "catch2/catch.hpp"
#include "utils/interp_2d.hh"
#include "utils/make_ifstream.hh"
#include "utils/read_vector.hh"

// This is the code we're actually testing
#include "models/nfw_sigma_mis.cuh"

#include <fstream>
#include <iomanip>
#include <vector>

using y3_cuda::NFW_SIGMA_MIS;
using quad::Interp2D;

// This test case (see https://github.com/catchorg/Catch2 for an introduction
// to the testing framework we are using for C++).
//
// This test case is not typical for Catch2 because it serves two purposes:
//   1. It is an actual unit test; see the use of the CHECK macro in the code
//      below.
//   2. As a side effect to running the test, a file is written that can be
//      used in creating a section for our validation document.
//


TEST_CASE("NFW profile function works")
{
  NFW_SIGMA_MIS model;
  auto out = y3_cluster::make_ofstream("data/nfw_sigma_mis_test.out");
  REQUIRE(out.good());
  out << std::setw(16);
  out << std::setprecision(16);
  out << "r\trmis\tlnM\tytrue\tytest\n";

  // Here we need to fill the vectors rs, rmiss, lnMs and answers with the
  // right data for the test.
  // CURRENTLY THIS CODE DOES NOT FILL THE ARRAYS WITH RIGHT DATA.
  std::vector<double> rs{1.0};
  std::vector<double> rmiss{0.5};
  std::vector<double> lnMs{10.0};
  std::vector<double> answers{1.0};
    
  double const epsrel = 1.0e-3;
  double const epsabs = 1.0e-12;

  std::size_t const n_rs = rs.size();
  std::size_t const n_rmiss = rmiss.size();
  std::size_t const n_lnMs = lnMs.size();

  for (std::size_t i = 0; i != n_rs; ++i) {
    for (std::size_t j = 0; j != n_rmiss; ++j) {
      for (std::size_t k = 0; k != n_lnMs; ++k) {
        double const calculated_result = model(rs[i], rmiss[j], lnMs[k]);
        // Adjust this look to look up the answers correctly; the indexing needed will depend
        // on the order of storage of the results.
        double const expected_result = answers[k + n_lnMs*j + n_lnMs*n_rmiss*i];
        CHECK(calculated_result == Approx(expected_result).epsilon(epsrel).margin(epsabs));
        out << rs[i] << '\t'
            << rmiss[j] << '\t'
            << lnMs[k] << '\t'
            << expected_result << '\t'
            << calculated_result << '\n';
      }
    }
  }
}
