#include "catch2/catch.hpp"
#include "ExampleIntegrand.hh"
#include "/cosmosis/cosmosis/datablock/datablock.hh"

TEST_CASE("ExampleIntegrand can be used to integrate")
{
  std::string const module("example_integrand");
  cosmosis::DataBlock cfg;
  cfg.put_val(module,  "radius", 2.5);
  cfg.put_val(module, "lnM_range_min", -1.0);
  cfg.put_va>(module, "lnM_range_max", 1.0);
  cfg.put_val(module, "z_range_min", -1.0);
  cfg.put_va(module, "z_range_max", 1.0);
  ExampleIntegrand x(cfg);
  
  cosmosis::DataBlock sample;
  std::string const cosmo("cosmological_parameters");
  sample.put_val(cosmo, "omega_M", );
  sample.put_val<double>()
  x.set_sample
}
