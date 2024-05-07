#include "catch2/catch.hpp"
#include "models/hmf_t.hh"
#include "utils/make_ifstream.hh"
#include "utils/read_vector.hh"

// This is the code we're actually testing
#include "models/nfw_sigma_mis.cuh"
#include "utils/cuda_interp_2d.cuh"

#include <fstream>
#include <iomanip>
#include <fmt/format.h>

using y3_cluster::HMF_t;
// using y3_cluster::Interp2D;
using y3_cuda::NFW_SIGMA_MIS;

TEST_CASE("Test NFW Misc Implementation")
{
  // Setup the initial conditions properly
  double  conc = 4.0;
  double  omega_m = 0.3;
  double  RHOC = 2.77533742639e+11;
  double  RHOM = RHOC*omega_m;
  NFW_SIGMA_MIS model(conc, RHOM, "gamma");

  double const epsrel = 1.0e-3;
  double const epsabs = 1.0e-12;

  // Here we need to fill the vectors rs, rmiss, lnMs and answers 
  // Radii vector log spaced from 0.1 to 10 with 10 points
  std::vector<double> r{ 0.1       ,  0.16681005,  0.27825594,  0.46415888,  0.77426368,
                         1.29154967,  2.15443469,  3.59381366,  5.9948425};

    // Answers drawn from cluster toolkit
  // Each columns is a different rmis, each row is a different r
  std::vector<double> answers{154.2477729, 36.2597318, 5.9036604,
                              126.1577857, 35.0800970, 5.8699027,
                               86.6157755, 32.5206407, 5.7884943,
                               47.3583660, 27.6699473, 5.6047940,
                               21.4019571, 20.1449034, 5.2255105,
                                8.8737342, 11.4269575, 4.5301542,
                                3.5299064, 4.7155369, 3.4469777,
                                1.3583749, 1.5653792, 2.1094762,
                                0.5101472, 0.5325816, 0.9306969,
                                0.1884599, 0.1911689};

  std::vector<double> rmis{0.1,0.5,2.0};
  std::vector<double> lnMs{32.2361913};

  std::size_t const n_r = r.size();
  std::size_t const n_rmiss = rmis.size();
  std::size_t const n_lnMs = lnMs.size();

  auto out = y3_cluster::make_ofstream("data/nfw_sigma_mis_test.out");
  REQUIRE(out.good());
  out << std::setw(16);
  out << std::setprecision(16);
  out << "r\trmis\tlnM\tytrue\tytest\n";
  
  for (std::size_t i = 0; i != n_r; ++i) {
    for (std::size_t j = 0; j != n_rmiss; ++j) {
      for (std::size_t k = 0; k != n_lnMs; ++k) {
        double const calculated_result = model(r[i], rmis[j], lnMs[k]);
        // Adjust this look to look up the answers correctly; the indexing needed will depend
        // on the order of storage of the results.
        double const expected_result = answers[k + n_lnMs*j + n_lnMs*n_rmiss*i];
        CHECK(calculated_result == Approx(expected_result).epsilon(epsrel).margin(epsabs));
        out << r[i] << '\t'
            << rmis[j] << '\t'
            << lnMs[k] << '\t'
            << expected_result << '\t'
            << calculated_result << '\n';
      }
    }
  }
}

// 
// To Generate The Answers; Go To a Python Interface with Cluster Toolkit
// 
// import numpy as np
// import cluster_toolkit as ct

// def make_sigma_mis(R, Rmis=0.3, Mt=1e14, c=4., omega_m=0.3):
//     Rp = np.logspace(np.log10(0.001),np.log10(15.),300)
//     st = ct.deltasigma.Sigma_nfw_at_R(Rp, Mt, c, omega_m)
//     _st_mis = ct.miscentering.Sigma_mis_at_R(Rp, Rp, st, Mt, c, omega_m, Rmis, kernel='gamma')
//     return np.interp(R,Rp,_st_mis)

// R = np.logspace(np.log10(0.1),np.log10(10.),10)
// R

// res = np.vstack([make_sigma_mis(R, 0.1),
//                  make_sigma_mis(R, 0.5),
//                  make_sigma_mis(R, 2.)])

// for row in res.T:
//     print('%.7f, %.7f, %.7f,'%(row[0],row[1],row[2]))