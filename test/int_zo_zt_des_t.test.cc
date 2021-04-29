#include "catch2/catch.hpp"
#include "models/int_zo_zt_des_t.hh"
#include "models/sigma_photoz_des.hh"
#include "utils/make_ifstream.hh"

#include <fstream>
#include <iostream>
#include <string>

using y3_cluster::INT_ZO_ZT_DES_t;
using y3_cluster::SIGMA_PHOTOZ_DES_t;

TEST_CASE("int_zo_zt_des and sigma_photoz_des objects works")
{
  std::ifstream infile =
    y3_cluster::make_ifstream("data/test_int_zo_zt_des.txt");

  // Test that we can open the file
  REQUIRE(infile.good());

  // Initialize vectors to hold truth
  std::vector<double> zmin, zmax, zt;
  std::vector<double> sigma_true, int_true;
  std::string headerline;

  // Read in the truth
  std::getline(infile, headerline);
  std::getline(infile, headerline);
  while (infile) {
    double onezmin, onezmax, onezt;
    double onesig, onepdf, oneint;
    infile >> onezmin >> onezmax >> onezt;
    infile >> onesig >> oneint;
    zmin.push_back(onezmin);
    zmax.push_back(onezmax);
    zt.push_back(onezt);
    sigma_true.push_back(onesig);
    int_true.push_back(oneint);
  }

  // Ensure that we have the some number of values in each column
  REQUIRE(zmin.size() == zmax.size());
  REQUIRE(zmax.size() == zt.size());
  REQUIRE(zt.size() == sigma_true.size());
  REQUIRE(zt.size() == int_true.size());

  // Create an instance of each class
  SIGMA_PHOTOZ_DES_t sigma_photoz_des;
  INT_ZO_ZT_DES_t intzozt;

  // Test the model against truth
  double constexpr epsrel = 1.0e-6;
  double constexpr epsabs = 1.0e-12;
  for (std::size_t i = 0, len = int_true.size(); i != len; ++i) {
    double const sigma_test = sigma_photoz_des(zt[i]);
    double const int_test = intzozt(zmin[i], zmax[i], zt[i]);
    CHECK(sigma_test == Approx(sigma_true[i]).epsilon(epsrel).margin(epsabs));
    CHECK(int_test == Approx(int_true[i]).epsilon(epsrel).margin(epsabs));
  }
}
