#include "catch2/catch.hpp"
#include "utils/multi_dimensional_integrator.hh"

#include <cstddef>
#include <stdexcept>


using y3_cluster::MultiDimensionalIntegrator;

TEST_CASE("can use all algorithms to integrate a constant function over the unit square")
{
  cubacores(0, 0);
  auto flat_func = [](double, double) { return 1.0; };

  MultiDimensionalIntegrator algs;
  double const epsrel = 1.0e-6;
  double const epsabs = 1.0e-12;
  for (std::size_t i = 0; i !=  algs.num_algs(); ++i)
  {
    auto res = algs.integrate(i, flat_func, epsrel, epsabs);
    CHECK(res.status == 0);
    CHECK(res.value == Approx(1.0).epsilon(epsrel));
  }
}

TEST_CASE("can use all algorithms to integrate a constant function over a specified area")
{
  cubacores(0, 0);
  auto flat_func = [](double, double) { return 1.0; };
  cubacpp::IntegrationVolume<2> vol {{-1.0, -1.0}, {  1.0, 1.0}};

  MultiDimensionalIntegrator algs;
  double const epsrel = 1.0e-6;
  double const epsabs = 1.0e-12;
  for (std::size_t i = 0; i !=  algs.num_algs(); ++i)
  {
    auto res = algs.integrate(i, flat_func, epsrel, epsabs, vol);
    CHECK(res.status == 0);
    CHECK(res.value == Approx(4.0).epsilon(epsrel));
  }
}

TEST_CASE("can use named algorithms to integrate a constant function over the unit square")
{
  cubacores(0, 0);
  auto flat_func = [](double, double) { return 1.0; };
  double const epsrel = 1.0e-6;
  double const epsabs = 1.0e-12;
  
  SECTION("using cuhre")
  {
    MultiDimensionalIntegrator algs("cuhre");
    auto res = algs.integrate(flat_func, epsrel, epsabs);
    CHECK(res.status == 0);
    CHECK(res.value == Approx(1.0).epsilon(epsrel));
  }

  SECTION("using vegas")
  {
    MultiDimensionalIntegrator algs("vegas");
    auto res = algs.integrate(flat_func, epsrel, epsabs);
    CHECK(res.status == 0);
    CHECK(res.value == Approx(1.0).epsilon(epsrel));
  }

  SECTION("using suave")
  {
    MultiDimensionalIntegrator algs("suave");
    auto res = algs.integrate(flat_func, epsrel, epsabs);
    CHECK(res.status == 0);
    CHECK(res.value == Approx(1.0).epsilon(epsrel));
  }
}

TEST_CASE("wrong algorithm name is illegal")
{
  REQUIRE_THROWS_AS(MultiDimensionalIntegrator("unknown"), std::runtime_error);
}