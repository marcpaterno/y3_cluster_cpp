#include "catch2/catch.hpp"
#include "a_cen_t.hh"
#include "t_cen_t.hh"

using y3_cluster::A_CEN_t;
using y3_cluster::T_CEN_t;

TEST_CASE("a_cen_t and t_cen_t works")
{
  constexpr std::size_t ntests = 12;
  // Unit test with the files: DS_M_2e+14_4e+14_z_0_0.34.dat and DS_M_2e+14_4e+14_z_0_0.34_cosi_0.8_1.dat
  // rp: Units of [Mpc/h]
  constexpr double rp[ntests] = {
      0.125893, 0.199526, 0.316228, 0.501187,
      0.794328, 1.25893, 1.99526, 3.16228,
      5.01187, 7.94328, 12.5893, 19.9526
  };
  // units of [h Msun/pc^2]
  constexpr double ds[ntests] = {
      167.484, 141.745, 109.689, 79.222,
      53.9142, 34.3401, 19.7749, 10.0391,
      4.9866, 2.79118, 1.71323, 1.09489
  };
  // units of [h Msun/pc^2]
  constexpr double dsmu[ntests] = {
      215.694, 179.938, 136.966, 97.6729,
      65.1235, 39.7234, 22.1475, 11.6005,
      6.06163, 3.45574, 2.09941, 1.32538
  };

  A_CEN_t acen;
  T_CEN_t tcen;

  for (std::size_t i = 0, sz = ntests; i != sz; ++i)
  {
    // For the bin 0.8<mu<1
    const double mu = 0.9;
    // TODO: The 0.0 is filling in for ln(Mass) right now.
    // Should there be a mass dependence? Zhang, Wu, and Zhang does not
    // list one, so if we are using their model we don't need it.
    double const fz = std::exp(acen(mu) * tcen(rp[i], 0.0));
    double constexpr epsrel = 1.0e-1;
    CHECK(fz == Approx(dsmu[i]/ds[i]).epsilon(epsrel));
  }
}
