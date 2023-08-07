#include "catch2/catch.hpp"
#include "models/hmf_t.hh"
#include "utils/make_ifstream.hh"
#include "utils/read_vector.hh"

// This is the code we're actually testing
#include "models/nfw_dsigma_mis.cuh"
#include "utils/cuda_interp_2d.cuh"

#include <fstream>
#include <iomanip>
#include <fmt/format.h>

using y3_cluster::HMF_t;
// using y3_cluster::Interp2D;
using y3_cuda::NFW_DSIGMA_MIS;

TEST_CASE("Test NFW Misc Implementation")
{
  // Setup the initial conditions properly
  double  conc = 4.0;
  double  omega_m = 0.3;
  double  RHOC = 2.77533742639e+11;
  double  RHOM = RHOC*omega_m*1e-12;
  NFW_DSIGMA_MIS model(conc, RHOM, "gamma");

  double const epsrel = 1.0e-3;
  double const epsabs = 1.0e-12;

  // Here we need to fill the vectors rs, rmiss, lnMs and answers with the
  // right data for the test.
  // CURRENTLY THIS CODE DOES NOT FILL THE ARRAYS WITH RIGHT DATA.
  std::vector<double> r{0.1       ,  0.10476158,  0.10974988,  0.1149757 ,  0.12045035,
                        0.12618569,  0.13219411,  0.13848864,  0.14508288,  0.15199111,
                        0.15922828,  0.16681005,  0.17475284,  0.18307383,  0.19179103,
                        0.2009233 ,  0.21049041,  0.22051307,  0.23101297,  0.24201283,
                        0.25353645,  0.26560878,  0.27825594,  0.29150531,  0.30538555,
                        0.31992671,  0.33516027,  0.35111917,  0.36783798,  0.38535286,
                        0.40370173,  0.42292429,  0.44306215,  0.46415888,  0.48626016,
                        0.5094138 ,  0.53366992,  0.55908102,  0.58570208,  0.61359073,
                        0.64280731,  0.67341507,  0.70548023,  0.7390722 ,  0.77426368,
                        0.81113083,  0.84975344,  0.89021509,  0.93260335,  0.97700996,
                        1.02353102,  1.07226722,  1.12332403,  1.17681195,  1.23284674,
                        1.29154967,  1.35304777,  1.41747416,  1.48496826,  1.55567614,
                        1.62975083,  1.70735265,  1.78864953,  1.87381742,  1.96304065,
                        2.05651231,  2.15443469,  2.25701972,  2.36448941,  2.47707636,
                        2.59502421,  2.71858824,  2.84803587,  2.98364724,  3.12571585,
                        3.27454916,  3.43046929,  3.59381366,  3.76493581,  3.94420606,
                        4.1320124 ,  4.32876128,  4.53487851,  4.75081016,  4.97702356,
                        5.21400829,  5.46227722,  5.72236766,  5.9948425 ,  6.28029144,
                        6.57933225,  6.8926121 ,  7.22080902,  7.56463328,  7.92482898,
                        8.30217568,  8.69749003,  9.11162756,  9.54548457, 10.        };

  std::vector<double> answers{0.13471281, 0.13461068, 0.13557133, 0.13754422, 0.14053754,
                              0.14454731, 0.14958792, 0.1562642 , 0.16275065, 0.171022  ,
                              0.18034949, 0.19087506, 0.20258338, 0.21554443, 0.23076525,
                              0.24529831, 0.26233365, 0.28079197, 0.30079245, 0.32244284,
                              0.34568505, 0.37079382, 0.39763714, 0.42629263, 0.45705754,
                              0.48993903, 0.52485361, 0.56193001, 0.60131041, 0.64309236,
                              0.68735942, 0.73396995, 0.78314185, 0.83490599, 0.88925749,
                              0.94625767, 1.00594054, 1.06827617, 1.13324641, 1.20080647,
                              1.27090306, 1.34422162, 1.41826402, 1.49558622, 1.57460229,
                              1.65588622, 1.73869635, 1.82286253, 1.90846995, 1.99483428,
                              2.08172062, 2.16869556, 2.25540928, 2.34159988, 2.42650972,
                              2.50975788, 2.59085978, 2.66920964, 2.74420511, 2.81544493,
                              2.8822743 , 2.94402435, 3.00021814, 3.05027111, 3.09365498,
                              3.12989022, 3.15851301, 3.17911654, 3.19136299, 3.19505296,
                              3.18990406, 3.17581201, 3.15276684, 3.12083323, 3.08021542,
                              3.03112707, 2.97395734, 2.90915796, 2.83718823, 2.75895686,
                              2.67483054, 2.5855551 , 2.49211001, 2.39523267, 2.29569984,
                              2.19437233, 2.09211888, 1.98972258, 1.88789133, 1.78740156,
                              1.68890786, 1.59299831, 1.50041633, 1.41165886, 1.3272736 ,
                              1.24787107, 1.17424161, 1.10705334, 1.04704371, 0.99389509};

  std::vector<double> rmis{1.0};
  std::vector<double> lnMs{14.0};

  std::size_t const n_r = r.size();
  std::size_t const n_rmiss = rmis.size();
  std::size_t const n_lnMs = lnMs.size();

  auto out = y3_cluster::make_ofstream("data/nfw_dsigma_mis_test.out");
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
// omega_m = 0.3 
// Mt = 1e14
// c = 4.
// Rmis = 1.0

// R = np.logspace(np.log10(0.1),np.log10(10.),100)
// st = ct.deltasigma.Sigma_nfw_at_R(R, Mt, c, omega_m)
// st_mis = ct.miscentering.Sigma_mis_at_R(R, R, st, Mt, c, omega_m, Rmis, kernel='gamma')
// dst_mis = ct.miscentering.DeltaSigma_mis_at_R(R, R, st_mis)