#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "modules/ExampleVectorIntegrand.hh"
#include "catch2/catch.hpp"
#include "cubacpp/cuhre.hh"

#include <iostream>
#include <vector>

TEST_CASE("2D example integrand")
{
  std::string const module_label("example_vector_integrand");
  cosmosis::DataBlock cfg;
  std::vector<double> radii{2.5, 5.0};
  cfg.put_val(module_label, "radii", radii);
  ExampleVectorIntegrand f(cfg);
  cosmosis::DataBlock sample;
  std::string const cosmo("cosmological_parameters");
  sample.put_val(cosmo, "sigma_8", 0.75);
  f.set_sample(sample);

  SECTION("evaluation")
  {
    std::vector<double> expected{3.0, 2.0};
    std::vector<double> fxy = f(5.0, 1.75);
    CHECK(fxy.size() == 2);
    for (std::size_t i = 0; i != fxy.size(); ++i) {
      CHECK(fxy[i] == Approx(expected[i]).epsilon(1.0e-12));
    }
  } // evaluation
  
  SECTION("integrate unit volume")
  {
    
    double const epsrel = 1.0e-12;
    double const epsabs = 1.0e-12;
    cubacpp::Cuhre alg;
    using iv_t = cubacpp::IntegrationVolume<2>;
    iv_t unitvolume {};
    auto res = alg.integrate(f, epsrel, epsabs, unitvolume);
    CHECK(res.converged());
    std::vector<double> expected { 0.345833, 0.245833 };
    for (std::size_t i = 0, sz = res.value.size(); i != sz; ++i) {
      CHECK(res.value[i] == Approx(expected[i]).epsilon(1.0e-3));
    }
  }

  SECTION("integrate non-default volumes")
  {
    double const epsrel = 1.0e-12;
    double const epsabs = 1.0e-9;
    cubacpp::Cuhre alg;
    cosmosis::DataBlock cfg2;
    cfg2.put_val(module_label, "x_low", std::vector<double>{-3.0});
    cfg2.put_val(module_label, "x_high", std::vector<double>{2.0});
    cfg2.put_val(module_label, "y_low", std::vector<double>{4.0});
    cfg2.put_val(module_label, "y_high", std::vector<double>{7.0});
    using iv_t = cubacpp::IntegrationVolume<2>;
    std::vector<iv_t> vols = ExampleVectorIntegrand::make_integration_volumes(cfg2);
    CHECK(vols.size() == 1);
    CHECK(vols[0].jacobian() == 15.0);
    auto res = alg.integrate(f, epsrel, epsabs, vols[0]);
    CHECK(res.converged());
    std::vector<double> expected { 346.688, 348.187 };
    for (std::size_t i = 0, sz = res.value.size(); i != sz; ++i) {
     CHECK(res.value[i] == Approx(expected[i]).epsilon(1.0e-3));
    }
  }
}
