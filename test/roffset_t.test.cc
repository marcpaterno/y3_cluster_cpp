#include "catch2/catch.hpp"
#include "test/roffset_t.hh"

#include <fstream>

using y3_cluster::ROFFSET_t;

TEST_CASE("roffset_t works")
{
    std::ifstream infile {"test_Rmis.txt"};
    // Use REQUIRE for immediate failure if we can't open the file.
    REQUIRE(infile.good());
    std::vector<double> rmis;
    std::vector<double> prob;
    std::string line;
    
    // Remove the two lines of header
    getline(infile, line);
    getline(infile, line);
    
    while (infile)
    {
        double r, p;
        infile >> r >> p;
        rmis.push_back(r);
        prob.push_back(p);
    }
    
    infile.close();
    
    // Remove the extra copy of the last line
    rmis.pop_back();
    prob.pop_back();
    
    // If the file is well-formed, we have the same number of z-values as
    // y(z)-values.
    REQUIRE(rmis.size() == prob.size());
    
    const double tau = 0.150000;
    ROFFSET_t roffset(tau);
    
    for (std::size_t i = 0, sz = rmis.size(); i != sz; ++i)
    {
        double const fz = roffset(rmis[i]);
        double constexpr epsrel = 1.0e-6;
        CHECK(fz == Approx(prob[i]).epsilon(epsrel));
    }
}
