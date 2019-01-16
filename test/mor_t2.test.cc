#include "catch2/catch.hpp"
#include "mor_t2.hh"

#include <fstream>
#include <iostream>

using y3_cluster::MOR_t2;

TEST_CASE("mor_t2 works")
{
    std::ifstream infile {"../data/mor_tt_test.txt"};
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
        log10m.push_back(lm*std::log(10));
        lambda_true.push_back(lt);
        prob.push_back(p);
    }
    
    infile.close();
    log10m.pop_back();
    lambda_true.pop_back();
    prob.pop_back();

    // If the file is well-formed, we have the same number of z-values as
    // y(z)-values.
    REQUIRE(log10m.size() == lambda_true.size());
    REQUIRE(log10m.size() == prob.size());

    const double sigma_intr = 0.10;
    const double alpha = 0.65;

    MOR_t2 mor_t(pow(10,11.2), pow(10,12.42), alpha, sigma_intr);
   
    for (std::size_t i = 0, sz = log10m.size(); i != sz; ++i)
    {
        double const fz = mor_t(lambda_true[i], log10m[i], 1.0);
        double constexpr epsrel = 5.0e-5;
        double constexpr epsabs = 1.0e-10;
        CHECK( fz == Approx(prob[i]).epsilon(epsrel).margin(epsabs) );
    }
}
