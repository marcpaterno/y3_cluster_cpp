#include "catch2/catch.hpp"
#include "models/sigma_photoz_des.hh"
#include "models/zo_zt_des_t.hh"

#include <fstream>
#include <string>
#include <iostream>

using y3_cluster::ZO_ZT_DES_t;
using y3_cluster::SIGMA_PHOTOZ_DES_t;

TEST_CASE("P(zobs|ztrue) implementation zo_zt_des works") {
    std::ifstream infile {"../data/test_zo_zt_des.txt"};

    // Test that we can open the file
    REQUIRE(infile.good());

    // Initialize vectors to hold truth
    std::vector<double> zo, zt, sig_true, prob_true;;
    std::string headerline;

    // Read in the truth
    std::getline(infile, headerline);
    std::getline(infile, headerline);
    while(infile) {
        double onezo, onezt, onesig, oneprob;
        infile >> onezo >> onezt >> onesig >> oneprob;
        zo.push_back(onezo);
        zt.push_back(onezt);
        sig_true.push_back(onesig);
        prob_true.push_back(oneprob);
    }

    // Ensure that we have the some number of values in each column
    REQUIRE(zo.size() == zt.size());
    REQUIRE(zt.size() == prob_true.size());
    REQUIRE(zt.size() == sig_true.size());
 
    // Create an instance of the class
    SIGMA_PHOTOZ_DES_t sigma_photoz_des;
    ZO_ZT_DES_t pzozt;

    // Test the model against truth
    double constexpr epsrel = 1.0e-6;
    double constexpr epsabs = 1.0e-12;
    for(std::size_t i=0, len=prob_true.size(); i!=len; ++i) {
        double const sigma_test = sigma_photoz_des(zt[i]);
        double const prob_test = pzozt(zo[i], zt[i]);
        CHECK(sigma_test == Approx(sig_true[i]).epsilon(epsrel).margin(epsabs));
        CHECK(prob_test == Approx(prob_true[i]).epsilon(epsrel).margin(epsabs));
    }
}
