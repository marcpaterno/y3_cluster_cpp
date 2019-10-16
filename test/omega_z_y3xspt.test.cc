#include "catch2/catch.hpp"
#include "models/sptxdes/omega_z_y3xspt.hh"

#include <fstream>
#include <iomanip>
#include <string>

using y3_cluster::OMEGA_Z_Y3XSPT;

TEST_CASE("omega_z_y3xspt works")
{
  std::ifstream infile{"../data/test_omega_z_y3xspt.txt"};

  // Test that we can open the file
  REQUIRE(infile.good());

  // Initialize vectors to hold truth
  std::vector<double> zs;
  std::vector<double> area_true;
  std::string headerline;

  // Read from the truth file
  std::getline(infile, headerline);
  std::getline(infile, headerline);
  while (infile) {
    double z, area;
    infile >> z >> area;
    zs.push_back(z);
    area_true.push_back(area);
  }

  // Ensure that we have the some number of values in each column
  REQUIRE(zs.size() == area_true.size());

  // Create an instance of the struct
  OMEGA_Z_Y3XSPT omega;

  // Prepare the output file for validation
  std::ofstream out{"../data/omega_z_y3xspt.out"};
  out << std::setw(16);
  out << std::setprecision(16);
  out << "z\tAtrue\tAtest\n";

  // Test computation against truth
  double constexpr epsrel = 1.0e-6;
  double constexpr epsabs = 1.0e-12;
  for (std::size_t i = 0, len = zs.size(); i != len; ++i) {
    double const area_test = omega(zs[i]);
    CHECK(area_test == Approx(area_true[i]).epsilon(epsrel).margin(epsabs));
    out << zs[i] << '\t' << area_true[i] << '\t' << area_test << '\n';
  }
}
