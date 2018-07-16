#include "catch2/catch.hpp"

#include <array>
#include <tuple>
#include <utility>

#include <integration_range.hh>
#include <test/param_space_explorer.hh>

using y3_cluster::IntegrationRange,
      std::array,
      std::tuple,
      std::pair;
TEST_CASE("Param space explorer acts as expected")
{
    double a, b;
    IntegrationRange air{3.0, 11.0},
                     bir{-13.4, 19.3};

    ParamSpaceExplorer<2> pse2d({&a, &b},
                                {air, bir},
                                {4, 5});

    array<pair<double, double>, 4*5> expected_2d;

    for (auto i = 0u; i < 4; i++)
        for (auto j = 0u; j < 5; j++)
            expected_2d[5*i + j] = {air.transform(((double) i) / 3),
                                    bir.transform(((double) j) / 4)};

    auto idx = 0u;
    do {
        auto [aexp, bexp] = expected_2d[idx++];
        CHECK(a == Approx(aexp).margin(1e-5));
        CHECK(b == Approx(bexp).margin(1e-5));
    } while (pse2d.next());


    double c, d, e;
    IntegrationRange cir{-100, -1000},
                     dir{-10, 3},
                     eir{1032, 9000};

    ParamSpaceExplorer<3> pse3d({&c, &d, &e},
                                {cir, dir, eir},
                                {13, 5, 2});

    array<tuple<double, double, double>, 13*5*2> expected_3d;

    for (auto i = 0u; i < 13; i++)
        for (auto j = 0u; j < 5; j++)
            for (auto k = 0u; k < 2; k++)
                expected_3d[5*2*i + j*2 + k] = {cir.transform(((double) i) / 12),
                                                dir.transform(((double) j) / 4),
                                                eir.transform(((double) k) / 1)};

    idx = 0u;
    do {
        auto [cexp, dexp, eexp] = expected_3d[idx++];
        CHECK(c == Approx(cexp).margin(1e-5));
        CHECK(d == Approx(dexp).margin(1e-5));
        CHECK(e == Approx(eexp).margin(1e-5));
    } while (pse3d.next());
}
