#include "utils/module_macros.hh"
#include "utils/make_integration_volumes.hh"

#include "cosmosis/datablock/datablock.hh"
#include "cubacpp/integration_result.hh"
#include "cubacpp/integration_volume.hh"

#include "models/hmf_t.hh"
#include "models/dv_do_dz_t.hh"
#include "models/omega_z_des.hh"
#include "models/y_sz.hh"

#include <iostream>
#include <optional>
#include <vector>

using namespace y3_cluster;

// compton_y_sims is a module to compute the expected Compton-y profile of
// galaxy clusters
//
// Notes:
//    1) std::optional<T> is used for data members that are not
//    constructible from CosmoSIS configuration parameters.
//
//    2) The object as created by the only constructor does not need to be
//    in a callable state.
//
//    3) After a call to set_sample has been made, the object must be in a
//    callable state.
//
//    4) State that *can* be correctly set by the constructor *should* be set by
//    the constructor. Do not needlessly repeat initialization in calls to
//    set_sample.
//
//
class compton_y_sims {
  // Integrate over M and z - 2d
  using volume_t = cubacpp::IntegrationVolume<2>;

private:
  // State obtained from each sample.
  // If there were a type X that did not have a default constructor,
  // we would use std::optional<X> as our data member.
  std::optional<OMEGA_Z_DES> omega_z;
  std::optional<DV_DO_DZ_t> dv_do_dz;
  std::optional<HMF_t> hmf;
  std::optional<Y_SZ> y_sz;
  // The angles (in arcmin) at which we want to compute expected Compton-y
  // profiles
  std::vector<double> thetas;
  // The low and high z bins
  // (needed for normalization)
  std::vector<double> z_low, z_high;

public:
  // Initialize my integrand object from the parameters read
  // from the relevant block in the CosmoSIS ini file.
  explicit compton_y_sims(cosmosis::DataBlock& config);

  // Set any data members from values read from the current sample.
  // Do not attempt to copy the sample!.
  void set_sample(cosmosis::DataBlock& sample);

  // The function to be integrated. All arguments to this function must be of
  // type double, and there must be at least two of them (because our
  // integration routine does not work for functions of one variable). The
  // function is const because calling it does not change the state of the
  // object.
  std::vector<double> operator()(double M, double z) const;

  // finalize_sample() is where products can be put into the cosmosis::DataBlock
  // representing the current sample. The object 'sample' passed to this
  // function will be the exact same object as was passed to the most recent
  // call to set_sample(). The object 'results' will be the results of the
  // integration that has just been done for that sample. This is generally the
  // item which should be put into the sample.
  void finalize_sample(
    cosmosis::DataBlock& sample,
    std::vector<cubacpp::integration_results_v> const& results) const;

  // module_label() should return the label for this module. The name this returns
  // is the name that must be used in the 'ini file' for configuring the module
  // made with this class.
  // We return char const* rather than std::string to avoid some needless memory
  // allocations.
  static char const* module_label();

  // The following non-member (static) function creates a vector of integration
  // volumes (the type alias defined above) based on the parameters read from
  // the configuration block for the module.
  static std::vector<volume_t> make_integration_volumes(
    cosmosis::DataBlock& cfg);
};

// We write using declarations so that we don't have to type the namespace name
// each time we use these names
using cosmosis::DataBlock;
using cubacpp::integration_results_v;

compton_y_sims::compton_y_sims(DataBlock& config)
    : omega_z()
    , dv_do_dz()
    , hmf()
    , y_sz()
    , thetas(get_datablock<std::vector<double>>(config, compton_y_sims::module_label(), "thetas"))
    , z_low(get_datablock<std::vector<double>>(config, compton_y_sims::module_label(), "z_low"))
    , z_high(get_datablock<std::vector<double>>(config, compton_y_sims::module_label(), "z_high"))
{
    if (z_low.size() != z_high.size())
        throw std::invalid_argument("compton_y_sims: z_low and z_high are different lengths");
}

void
compton_y_sims::set_sample(DataBlock& sample)
{
  // If we had a data member of type std::optional<X>, we would set the
  // value using std::optional::emplace(...) here. emplace takes a set
  // of arguments that it passes to the constructor of X.
  dv_do_dz.emplace(sample);
  hmf.emplace(sample);
  omega_z.emplace(sample);
  y_sz.emplace(sample);
}

std::vector<double>
compton_y_sims::operator()(double lnM, double z) const
{
  std::vector<double> results;
  results.reserve(thetas.size() + 1);

  // Number counts followed by the radius bins, repeating over zo bins
  double common_term = (*dv_do_dz)(z)
                       * (*hmf)(lnM, z)
                       * (*omega_z)(z);

  // This is the N integral
  results.push_back(common_term);

  for (auto ti = 0u; ti < thetas.size(); ti++) {
    // TODO make this real
    const double compton_y = (*y_sz)(std::exp(lnM), z, ti);
    results.push_back(common_term * compton_y);
  }

  return results;
}

void
compton_y_sims::finalize_sample(
  cosmosis::DataBlock& sample,
  std::vector<integration_results_v> const& results) const
{
  std::vector<int> NM_status;
  std::vector<int> NM_nregions;
  std::vector<int> NM_nevals;

  std::vector<double> N_vals;
  std::vector<double> N_errors;
  std::vector<double> N_probs;

  std::vector<double> compton_y_vals_temp;
  std::vector<double> compton_y_errors_temp;
  std::vector<double> compton_y_probs_temp;

  // TODO: Try to refactor this code, to abstract away the manipulations done for
  // each physics quantity.
  std::size_t n_radii_bins = thetas.size();
  const auto block_size = n_radii_bins + 1;

  if (z_low.size() != results.size())
    throw std::runtime_error("JOD - im wrong about results size");

  for (auto i = 0u; i < results.size(); i++) {
    auto const& result = results[i];

    if (result.value.size() != block_size)
      throw std::runtime_error("JOD - im wrong about results size");

    NM_status.push_back(result.status);
    NM_nregions.push_back(result.nregions);
    NM_nevals.push_back(result.neval);

    N_vals.push_back(result.value[0]);
    N_errors.push_back(result.error[0]);
    N_probs.push_back(result.prob[0]);

    // The output compton-y is weighted by the number counts
    // So we have to divide by the number counts
    for (auto ti = 1u; ti < result.value.size(); ++ti) {
      compton_y_vals_temp.push_back(result.value[ti] / result.value[0]);
      compton_y_errors_temp.push_back(result.error[ti] / result.value[0]);
      compton_y_probs_temp.push_back(result.prob[ti]);
    }
  }

  std::size_t expected_size = results.size() * n_radii_bins;
  if ((compton_y_vals_temp.size() != expected_size)
          || (compton_y_errors_temp.size() != expected_size)
          || (compton_y_probs_temp.size() != expected_size))
    throw std::runtime_error("internal error: dimension mismatch for output");

  std::vector<std::size_t> extents{results.size(), n_radii_bins};
  cosmosis::ndarray<double> compton_y_vals(compton_y_vals_temp, extents);
  cosmosis::ndarray<double> compton_y_errors(compton_y_errors_temp, extents);
  cosmosis::ndarray<double> compton_y_probs(compton_y_probs_temp, extents);

  sample.put_val(module_label(), "N_vals", N_vals);
  sample.put_val(module_label(), "N_errors", N_errors);
  sample.put_val(module_label(), "N_probs", N_probs);
  sample.put_val(module_label(), "NM_status", NM_status);
  sample.put_val(module_label(), "NM_nregions", NM_nregions);
  sample.put_val(module_label(), "NM_nevals", NM_nevals);

  sample.put_val(module_label(), "compton_y_radius_arcmin", thetas);

  sample.put_val(module_label(), "compton_y_vals", compton_y_vals);
  sample.put_val(module_label(), "compton_y_errors", compton_y_errors);
  sample.put_val(module_label(), "compton_y_probs", compton_y_probs);
}

char const*
compton_y_sims::module_label()
{
  return "compton_y_sims";
}

// The implementation of make_integration_volumes can be almost the same for
// any CosmoSISVectorIntegrand-type class. Only the names and number of the parameters
// provided need to be changed. It is critical that the names be given in the
// order that correspond to the order of arguments in the class's function call
// operator. While the compiler can verify the number of arguments provided is
// correct, it can not verify that their order matches the order of arguments in
// the function call operator.
std::vector<compton_y_sims::volume_t>
compton_y_sims::make_integration_volumes(cosmosis::DataBlock& cfg)
{
  try {
    const auto vols = y3_cluster::make_integration_volumes(cfg, module_label(), "M", "z");

    // We want to convert Ms to ln(Ms)
    std::vector<compton_y_sims::volume_t> corrected_vols;
    for (auto vol : vols) {
      auto [Mlo, zlo] = vol.transform({0.0, 0.0});
      auto [Mhi, zhi] = vol.transform({1.0, 1.0});

      compton_y_sims::volume_t corrected_vol{{std::log(Mlo), zlo},
                                             {std::log(Mhi), zhi}};
      corrected_vols.push_back(corrected_vol);
    }
    return corrected_vols;
  }
  catch (std::exception const& ex) {
    std::cerr << ex.what();
    throw;
  };
}

DEFINE_COSMOSIS_VECTOR_INTEGRATION_MODULE(compton_y_sims);
