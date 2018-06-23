#include "catch2/catch.hpp"
#include "cubacpp/cubacpp.hh"
#include "integration_range.hh"
#include "test/lc_lt_t.hh"

#include <fstream>
#include <iostream>
#include <string>

using y3_cluster::LC_LT_t, y3_cluster::IntegrationRange;
TEST_CASE("Lc_Lt_t works")
{
  std::ifstream infile{"test_lc_lt_t_SDSS.txt"};
  // Use REQUIRE for immediate failure if we can't open the file.
  REQUIRE(infile.good());
  std::vector<double> zs;
  std::vector<double> lts;
  std::vector<double> lcs;
  std::vector<double> ys;
  while (infile) {
    // We aren't bothering to test that the reading worked, because we're
    // careful to make sure the data file is not mal-formed when we write the
    // test.
    double z, lt, lc, y;
    infile >> z >> lt >> lc >> y;
    zs.push_back(z);
    lts.push_back(lt);
    lcs.push_back(lc);
    ys.push_back(y);
  }
  // If the file is well-formed, we have the same number of z-values as
  // y(z)-values.
  REQUIRE(zs.size() == lts.size());
  REQUIRE(lts.size() == lcs.size());
  REQUIRE(lcs.size() == ys.size());

  LC_LT_t lc_lt;

  for (std::size_t i = 0, sy = ys.size(); i != sy; ++i) {
    double const fz = lc_lt(lcs[i], lts[i], zs[i]);
    double constexpr epsrel = 1.0e-6;
    double constexpr epsabs = 1.0e-12;
    CHECK(fz == Approx(ys[i]).epsilon(epsrel).margin(epsabs));
  }

  // Lc_lt calculates $P(\lambda^{cen} | \lambda^{true}, z^{true})$ - so,
  // integral of Lc_lt over $\lambda^{cen}$ = [0, +inf] should be unity
  // Since total probability = 1
  // In practice - $\lambda^{cen}$ from 1 to 200 should be good enough
  IntegrationRange lc_ir{1.0, 200.0};
  IntegrationRange lt_ir{1.0, 100.0};
  IntegrationRange zt_ir{0.1, 0.3};

  cubacpp::Vegas v;
  v.maxeval = 9999999;
  const double epsrel = 1.0e-4;
  const double epsabs = 1.0e-12;

  const std::size_t width = 2;
  for (auto i = 0u; i < width; i++) {
    for (auto j = 0u; j < width; j++) {
      const double lt = lt_ir.transform((i + 1) / ((double)width + 1));
      const double zt = zt_ir.transform((j + 1) / ((double)width + 1));
      const auto res = v.integrate(
        [lt, zt, lc_ir, lc_lt](double scaled_lc) {
          const double lc = lc_ir.transform(scaled_lc);
          return lc_ir.jacobian() * lc_lt(lc, lt, zt);
        },
        epsrel,
        epsabs);

      // In reality - since we are integration [1, 200] not [0, +inf] - it
      // won't be exactly 1.0. So, arbitrary wiggle room 1e-2
      // TODO: should this epsilon be changed?
      CHECK(res.status == 0);
      CHECK(res.value == Approx(1.0).epsilon(2e-3).margin(2e-3));
    }
  }
}
