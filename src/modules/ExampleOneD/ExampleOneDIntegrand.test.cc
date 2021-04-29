#include "ExampleOneDIntegrand.hh"
#include "catch2/catch.hpp"
#include "cosmosis/datablock/datablock.hh"

#include <array>
#include <vector>

// Declare the funcs that the module exports
extern "C" {
void* setup(cosmosis::DataBlock*);
DATABLOCK_STATUS execute(cosmosis::DataBlock*, void*);
int cleanup(void*);
}

TEST_CASE("1D example integrand")
{
  SECTION("evaluation")
  {
    std::string const module_label("example_oned_integrand");
    cosmosis::DataBlock cfg;
    cfg.put_val(module_label, "b", 3.0);
    ExampleOneDIntegrand f(cfg);
    cosmosis::DataBlock sample;
    sample.put_val(module_label, "c", 1.2);
    f.set_sample(sample);
    f.set_grid_point({0.5});

    std::vector<double> xvals{-1.0, 0.0, 1.0, 3.33};
    std::vector<double> expected{3.7, 4.2, 4.7, 5.865};
    for (std::size_t i = 0; i != xvals.size(); ++i)
      CHECK(f(xvals[i]) == Approx(expected[i]).epsilon(1.0e-12));
  } // evaluation

  SECTION("test the integration")
  {
    std::string const module_label("example_oned_integrand");
    // Create the module.
    cosmosis::DataBlock cfg;
    cfg.put_val(module_label, "a", std::vector<double>{0.5, 2.0});
    cfg.put_val(module_label, "b", 3.0);
    cfg.put_val(module_label, "x_low", std::vector<double>{-3.0, -4.0});
    cfg.put_val(module_label, "x_high", std::vector<double>{2.0, 1.0});
    cfg.put_val(module_label, "eps_rel", 1.0e-10);
    cfg.put_val(module_label, "eps_abs", 1.0e-12);
    cfg.put_val(module_label, "max_eval", 250 * 1000);
    void* mod = setup(&cfg);
    CHECK(mod);

    // Execute the module once.
    cosmosis::DataBlock sample;
    sample.put_val(module_label, "c", 1.2);
    int rc = execute(&sample, mod);
    CHECK(rc == 0);

    // Check outputs that were saved to the datablock
    auto status = sample.view<std::vector<double>>(module_label, "status");
    auto vals = sample.view<std::vector<double>>(module_label, "vals");
    auto errors = sample.view<std::vector<double>>(module_label, "errors");
    std::vector<double> truth{19.75, 6.0};
    for (std::size_t i = 0; i != vals.size(); ++i) {
      CHECK(status[i] == 0);
      CHECK(vals[i] == Approx(truth[i]).epsilon(1.0e-12));
    }

    // Delete the module.
    rc = cleanup(mod);
    CHECK(rc == 0);
  }
}
