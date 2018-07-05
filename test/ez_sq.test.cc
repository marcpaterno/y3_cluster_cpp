#include "catch2/catch.hpp"
#include "test/ez_sq.hh"

#include <fstream>

using y3_cluster::EZ_sq;

double ez_sq_valid(double z, double omega_m, double omega_l, double omega_k)
{
    double const zplus1 = 1.0 + z;
    return (omega_m * zplus1 + omega_k ) * zplus1 * zplus1 +
              omega_l;
}

TEST_CASE("ez_sq works")
{
    double omega_m, omega_l, omega_k;
    omega_m = 0.27;
    omega_l = 0.63;
    omega_k = 1. - omega_m -omega_l;

    EZ_sq ez_sq(omega_m, omega_l, omega_k);

    double z = 0;
    const double dz = 0.01;
    const int nmax = 10000;

    double constexpr epsrel = 1.0e-6;

    for(int i = 0; i < nmax; i++)
    {
        double valid = ez_sq_valid(z, omega_m, omega_l, omega_k);
        CHECK(ez_sq(z) == Approx(valid).epsilon(epsrel));
        z += dz;
    }

}
