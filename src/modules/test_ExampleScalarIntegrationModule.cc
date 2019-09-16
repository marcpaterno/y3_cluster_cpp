#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"

// To test the module, we need to declare the functions
// that the module library exports. This is because modules
// have no headers.

extern "C" {
  void* setup(cosmosis::DataBlock*);
  DATABLOCK_STATUS execute(cosmosis::DataBlock*, void*);
  int cleanup(void*);
}

int main()
{
  std::string const module_label("example_scalar_integrand");
  // Create the module.
  cosmosis::DataBlock cfg;
  std::vector<double> radii{2.5, 5.0};
  cfg.put_val(module_label, "radii", radii);
  cfg.put_val(module_label, "x_low", std::vector<double>{-3.0});
  cfg.put_val(module_label, "x_high", std::vector<double>{2.0});
  cfg.put_val(module_label, "y_low", std::vector<double>{4.0});
  cfg.put_val(module_label, "y_high", std::vector<double>{7.0});
  cfg.put_val(module_label, "eps_rel", 1.0e-10);
  cfg.put_val(module_label, "eps_abs", 1.0e-12);
  cfg.put_val(module_label, "max_eval", 250*1000);  
  void* mod = setup(&cfg);
  if (!mod) return 1;

  // Execute the module once.
  cosmosis::DataBlock sample;
  std::string const cosmo("cosmological_parameters");
  sample.put_val(cosmo, "sigma_8", 0.75);
  int rc = execute(&sample, mod);
  if (rc != 0)  return 2;

  // Check values in 'sample' here.
  auto vals = sample.view<cosmosis::ndarray<double>>(module_label, "vals");

  // Delete the module.
  rc = cleanup(mod);
  if (rc != 0) return 4;
}

