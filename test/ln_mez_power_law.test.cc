#include "catch2/catch.hpp"
#include "models/ln_mez_power_law.hh"

#include <fstream>
#include <string>
#include <iomanip>
#include <iostream>

using y3_cluster::ln_mez_power_law;
TEST_CASE("ln_mez_power_law works")
{
    std::ifstream infile {"../data/test_ln_mez_power_law.txt"};

    // Test that we can open the file
    REQUIRE(infile.good());

    // Initialize vectors to hold truth
    std::size_t const ncols = 11;
    std::vector<std::vector<double> > inputs;
    std::vector<double> col(ncols, -1.0);
    std::string headerline;

    // Read in the truth
    std::getline(infile, headerline);
    std::getline(infile, headerline);
    while(infile) {
        for(std::size_t i=0; i!=ncols; ++i) {
            infile >> col[i];
        }
        inputs.push_back(col);
    }

    // Test code against truth
    double constexpr epsrel = 1.0e-6;
    double constexpr epsabs = 1.0e-12;
    for(std::size_t i=0, len=inputs.size(); i!=len; ++i) {
        double const A = inputs[i][0];
        double const B = inputs[i][1];
        double const C = inputs[i][2];
        double const lnM = inputs[i][3];
        double const lnMp = inputs[i][4];
        double const z = inputs[i][5];
        double const zp = inputs[i][6];
        double const omega_m = inputs[i][7];
        double const omega_l = inputs[i][8];
        double const omega_k = inputs[i][9];
        double const truth = inputs[i][10];
        ln_mez_power_law powlaw(A, B, C, lnMp, zp, omega_m, omega_l, omega_k);
        double const testval = powlaw(lnM, z);
        std::cout << truth << '\t' << testval << '\n';
        CHECK(testval == Approx(truth).epsilon(epsrel).margin(epsabs));
    }
}
