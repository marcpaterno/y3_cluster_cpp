#include "catch2/catch.hpp"

#include "sub_integral_memoizer.hh"

using y3_cluster::MemoizedIntegration;

TEST_CASE("Test integral memoization works") {
    MemoizedIntegration simple_polynomial([](double a, double x, double y) {
                const double a2 = 2*x + y,
                             a1 = 3*x*y - 3,
                             a0 = x/(y + 1) + x*x*4;
                return a2*a*a + a1*a + a0;
            },
            {0.0, 1.0},
            {0.0, 1.0}, 11,
            {0.0, 1.0}, 11);

    // Test for grid points
    for (auto i = 0; i < 11; i++) {
        for (auto j = 0; j < 11; j++) {
            const double x = ((double) i) / 10,
                         y = ((double) j) / 10,
                         a2 = 2*x + y,
                         a1 = 3*x*y - 3,
                         a0 = x/(y + 1) + x*x*4;

            CHECK(simple_polynomial(x, y) == Approx((1.0/3.0) * a2 + 0.5 * a1 + a0));
        }
    }

    for (auto i = 0; i < 10; i++) {
        for (auto j = 0; j < 10; j++) {
            const double x = ((double) 2*i + 1) / 20,
                         y = ((double) 2*j + 1) / 20,
                         a2 = 2*x + y,
                         a1 = 3*x*y - 3,
                         a0 = x/(y + 1) + x*x*4;

            CHECK(simple_polynomial(x, y) == Approx((1.0/3.0) * a2 + 0.5 * a1 + a0).margin(2e-2));
        }
    }
}
