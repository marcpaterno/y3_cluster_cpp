#include "catch2/catch.hpp"
#include "test/mor_t.hh"

#include "mz_power_law.hh"

#include <fstream>
#include <iostream>

using y3_cluster::MOR_t;
using y3_cluster::mz_power_law;

TEST_CASE("mor_t works")
{
    std::ifstream infile {"mor_tt_test.txt"};
    // Use REQUIRE for immediate failure if we can't open the file.
    REQUIRE(infile.good());
    std::vector<double> log10m;
    std::vector<double> lambda_true;
    std::vector<double> prob;
    std::string line;
    
    while (infile)
    {
        double lm, lt, p;
        infile >> lm >> lt >>p;
        log10m.push_back(lm);
        lambda_true.push_back(lt);
        prob.push_back(p);

    }
    
    infile.close();
    
    // Remove the extra copy of the last line
    log10m.pop_back();
    lambda_true.pop_back();
    prob.pop_back();

    // If the file is well-formed, we have the same number of z-values as
    // y(z)-values.
    REQUIRE(log10m.size() == lambda_true.size());
    REQUIRE(log10m.size() == prob.size());

    const double sigma = 0.10; // This is sigma_intr, not sigma. Fix it.
    const double alpha = 0.65;

    MOR_t mor_t(mz_power_law{1.e-14, 1., 0.1},sigma,alpha);
    
    for (std::size_t i = 0, sz = log10m.size(); i != sz; ++i)
    {
        double const fz = mor_t(lambda_true[i], log10m[i], 1.0);
        std::cout << prob[i] << " " << fz << std::endl;
        double constexpr epsrel = 1.0e-6;
        CHECK(fz == Approx(prob[i]).epsilon(epsrel));
    }
}
