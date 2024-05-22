#include "catch2/catch.hpp"
#include "utils/make_ifstream.hh"
#include "utils/read_vector.hh"

// This is the code we're actually testing
#include "models/int_lc_lt_des_t2.cuh"
#include "utils/cuda_interp_2d.cuh"

#include <fstream>
#include <iomanip>
#include <fmt/format.h>

// using y3_cluster::Interp2D;
using y3_cuda::INT_LC_LT_DES_t2;

TEST_CASE("Test Int Lc Lt DES t2 Implementation")
{
  // Setup the initial conditions properly
  INT_LC_LT_DES_t2 lc_lt;

  double const epsrel = 1.0e-5;
  double const epsabs = 1.0e-12;

  // Here we need to fill the vectors rs, rmiss, lnMs and answers 
  // Radii vector log spaced from 0.1 to 10 with 10 points
  std::vector<double> lc{21, 31, 41, 67};
  std::vector<double> zt{0.2,0.5,0.6};

  // Answers drawn from cluster toolkit
  // Each columns is a different rmis, each row is a different r
  std::vector<double> answers{0.58784, 0.548563, 0.518381, 0.65423, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

  std::size_t const n_lc = lc.size();
  std::size_t const n_zt = zt.size();

  auto out = y3_cluster::make_ofstream("data/int_lc_lt_des_t.out");
  REQUIRE(out.good());
  out << std::setw(16);
  out << std::setprecision(16);
  out << "r\trmis\tlnM\tytrue\tytest\n";
  
  for (std::size_t i = 0; i != n_lc; ++i) {
    for (std::size_t j = 0; j != n_zt; ++j) {
      double const calculated_result = lc_lt(lc[i], lc[i], zt[j]);
      // Adjust this look to look up the answers correctly; the indexing needed will depend
      // on the order of storage of the results.
      double const expected_result = answers[n_zt*j + n_lc*n_zt*i];
      CHECK(calculated_result == Approx(expected_result).epsilon(epsrel).margin(epsabs));
      out << lc[i] << '\t'
          << zt[j] << '\t'
          << expected_result << '\t'
          << calculated_result << '\n';
    }
  }
}