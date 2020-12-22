#include <cmath>
#include <iostream>
#include <memory>

#include "catch2/catch.hpp"
#include "cubacpp/gsl.hh"
#include "gsl/gsl_sf_bessel.h"
#include "models/dv_do_dz_t.hh"
#include "models/ez.hh"
#include "utils/integration_range.hh"
#include "utils/interp_1d.hh"
#include "utils/read_vector.hh"

TEST_CASE(
  "Check that different ways to integrate bessel functions are equivalent")
{
  double const h = 0.771358152;
  auto const zz = read_vector("distances/z.txt");
  auto const da_arr = read_vector("distances/d_a.txt");
  auto const dm_arr = read_vector("distances/d_m.txt");
  auto da_f = std::make_shared<y3_cluster::Interp1D const>(zz, da_arr);

  // Same cosmological parameters as used in the data above
  y3_cluster::EZ ez{0.1875, 0.81248, 0};
  y3_cluster::DV_DO_DZ_t dvdodz{da_f, ez, h};

  y3_cluster::Interp1D dcom(zz, dm_arr), inverse_dcom(dm_arr, zz);

  std::array<y3_cluster::IntegrationRange, 4> zirs{
    {{0.1, 0.3}, {0.3, 0.5}, {0.5, 0.7}, {0.7, 0.9}}};

  cubacpp::QAG qag(0, 0, GSL_INTEG_GAUSS61, 20);
  for (auto const& zir : zirs) {
    double const zstart = zir.transform(0.0), zend = zir.transform(1.0),
                 rstart = dcom.eval(zstart) * h, rend = dcom.eval(zend) * h;

    auto const r_normalization =
      (1.0 / 3.0) * (rend * rend * rend - rstart * rstart * rstart);

    auto const z_normalization =
      qag.with_range(zstart, zend)
        .integrate([&](double z) { return dvdodz(z); }, 1e-5, 1e-18);
    REQUIRE(z_normalization.status == 0);
    REQUIRE(r_normalization == Approx(z_normalization.value).epsilon(2e-4));

    y3_cluster::IntegrationRange ln_kir(std::log(0.0001), std::log(0.01));
    for (auto l = 0u; l < 10; l++) {
      for (auto ki = 0u; ki < 10; ki++) {
        double const lnk = ln_kir.transform(ki / 9.0), k = std::exp(lnk);

        auto const z_bessel_integral =
          qag.with_range(zstart, zend)
            .integrate(
              [&](double z) {
                return dvdodz(z) * gsl_sf_bessel_jl(l, k * dcom.eval(z) * h);
              },
              1e-5,
              1e-18);
        auto const r_bessel_integral =
          qag.with_range(rstart, rend)
            .integrate(
              [&](double r) { return r * r * gsl_sf_bessel_jl(l, k * r); },
              1e-5,
              1e-18);

        REQUIRE(z_bessel_integral.status == 0);
        REQUIRE(r_bessel_integral.status == 0);

        double const z_value = z_bessel_integral.value / z_normalization.value,
                     r_value = r_bessel_integral.value / r_normalization;
        CHECK(z_bessel_integral.value ==
              Approx(r_bessel_integral.value).epsilon(2e-3));
        CHECK(z_value == Approx(r_value).epsilon(1e-3).margin(1e-6));
      }
    }
  }
}
