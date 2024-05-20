#include "catch2/catch.hpp"
#include "models/mor_sdss_t.hh"
#include "utils/make_ifstream.hh"

#include <fstream>
#include <iomanip>

using y3_cluster::MOR_sdss;

TEST_CASE("mor_sdss works")
{
  std::ifstream infile = y3_cluster::make_ifstream("data/mor_tt_test.txt");
  // Use REQUIRE for immediate failure if we can't open the file.
  REQUIRE(infile.good());
  std::vector<double> log10m;
  std::vector<double> lambda_true;
  std::vector<double> prob;
  std::string line;

  while (infile) {
    double lm, lt, p;
    infile >> lm >> lt >> p;
    log10m.push_back(lm * std::log(10));
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

  MOR_sdss mor_t(pow(10, 11.2), pow(10, 12.42), alpha, sigma_intr);

  std::ofstream out = y3_cluster::make_ofstream("validation-data/mor_tt_test.out");
  REQUIRE(out.good());
  out << std::setw(16);
  out << std::setprecision(16);
  out << "lambda\tlog10m\tprobtrue\tprobtest\n";
  for (std::size_t i = 0, sz = log10m.size(); i != sz; ++i) {
    double const fz = mor_t(lambda_true[i], log10m[i], 1.0);
    double constexpr epsrel = 5.0e-6;
    double constexpr epsabs = 1.0e-10;
    CHECK(fz == Approx(prob[i]).epsilon(epsrel).margin(epsabs));
    out << lambda_true[i] << '\t' << log10m[i] << '\t' << prob[i] << '\t' << fz
        << '\n';
  }
}
