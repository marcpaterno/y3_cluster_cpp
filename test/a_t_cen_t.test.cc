#include <vector>

#include "catch2/catch.hpp"
#include "read_vector.hh"
#include "a_cen_t.hh"
#include "t_cen_t.hh"

using y3_cluster::A_CEN_t;
using y3_cluster::T_CEN_t;
using std::string;
using std::vector;
using std::pair;
using doubles = vector<double>;


// Reads a file of radius - DeltaSigma values, and returns two vectors:
//  { <radii>, <DeltaSigmas> }
pair<doubles, doubles>
read_DS_file(const string& filename)
{
  auto const data = read_vector(filename.c_str());
  CHECK((data.size() % 2) == 0);
  doubles radii, ds;
  for (auto i = 0u; i < data.size(); i += 2) {
    radii.push_back(data[i]);
    ds.push_back(data[i + 1]);
  }
  return { radii, ds };
}

TEST_CASE("a_cen_t and t_cen_t match simulation data")
{
  const vector<pair<double, string>> mass_bins{{
          {7.5e13, "DeltaSigma_for_Maria/DS_M_5e+13_1e+14_z_0_0.34"},
          {1.5e14, "DeltaSigma_for_Maria/DS_M_1e+14_2e+14_z_0_0.34"},
          {3e14, "DeltaSigma_for_Maria/DS_M_2e+14_4e+14_z_0_0.34"},
          {5.2e15, "DeltaSigma_for_Maria/DS_M_4e+14_1e+16_z_0_0.34"}}};
  const vector<pair<double, string>> mu_bins{{
          {0.1, "_cosi_0_0.2.dat"},
          {0.3, "_cosi_0.2_0.4.dat"},
          {0.5, "_cosi_0.4_0.6.dat"},
          {0.6, "_cosi_0.6_0.8.dat"},
          {0.9, "_cosi_0.8_1.dat"}}};

  A_CEN_t acen;
  T_CEN_t tcen;

  for (const auto [mass, mbinname] : mass_bins) {
    const auto [rs, ds] = read_DS_file(mbinname + ".dat");
    for (const auto [mu_mid, tail] : mu_bins) {
      const auto [rsmu, dsmu] = read_DS_file(mbinname + tail);
      CHECK(rs.size() == rsmu.size());

      for (auto i = 0u; i < rs.size(); ++i) {
        CHECK(rs[i] == Approx(rsmu[i]).epsilon(1e-12));
        double const fz = std::exp(acen(mu_mid) * tcen(rs[i], std::log(mass)));
        // TODO: This precision causes 4 assertions to fail!
        double constexpr epsrel = 1.0e-1;
        CHECK(fz == Approx(dsmu[i]/ds[i]).epsilon(epsrel));
      }
    }
  }
}
