#include "utils/make_integration_volumes.hh"
#include "utils/module_macros.hh"

#include "cosmosis/datablock/datablock.hh"
#include "cubacpp/integration_result.hh"
#include "cubacpp/integration_volume.hh"

#include "models/dv_do_dz_t.hh"
#include "models/hmf_t.hh"
#include "models/omega_z_sdss.hh"
#include "models/sig_sum.hh"

#include <optional>
#include <vector>
using namespace y3_cluster;

// buzzard_sigma_halos is a class that models the concept of
// "CosmoSISIntegrand", and is thus suitable for use as the template parameter
// for the class template CosmosisIntegrationModule.
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
class buzzard_sigma_halos {
private:
  // We define the type alias volume_t to be the right dimensionality
  // of integration volume for our integrand. If we were to change the
  // number of arguments required by the function call operator (below),
  // we would need to also modify this type alias to keep consistent.
  using volume_t = cubacpp::IntegrationVolume<2>;

  // State obtained from each sample.
  // If there were a type X that did not have a default constructor,
  // we would use std::optional<X> as our data member.
  std::optional<OMEGA_Z_SDSS> omega_z_sdss;
  std::optional<DV_DO_DZ_t> dv_do_dz;
  std::optional<HMF_t> hmf;
  std::optional<SIG_SUM> sigma;
  std::vector<double> radii_;

public:
  // Initialize my integrand object from the parameters read
  // from the relevant block in the CosmoSIS ini file.
  explicit buzzard_sigma_halos(cosmosis::DataBlock& config);

  // Set any data members from values read from the current sample.
  // Do not attempt to copy the sample!.
  void set_sample(cosmosis::DataBlock& sample);

  // The function to be integrated. All arguments to this function must be of
  // type double, and there must be at least two of them (because our
  // integration routine does not work for functions of one variable). The
  // function is const because calling it does not change the state of the
  // object.
  std::vector<double> operator()(double zt, double lnM) const;

  // finalize_sample() is where products can be put into the cosmosis::DataBlock
  // representing the current sample. The object 'sample' passed to this
  // function will be the exact same object as was passed to the most recent
  // call to set_sample(). The object 'results' will be the results of the
  // integration that has just been done for that sample. This is generally the
  // item which should be put into the sample.
  void finalize_sample(
    cosmosis::DataBlock& sample,
    std::vector<cubacpp::integration_results_v> const& results) const;

  // module_label() should return the label for this module. The name this
  // returns is the name that must be used in the 'ini file' for configuring the
  // module made with this class. We return char const* rather than std::string
  // to avoid some needless memory allocations.
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

buzzard_sigma_halos::buzzard_sigma_halos(DataBlock& config)
  : omega_z_sdss()
  , dv_do_dz()
  , hmf()
  , sigma()
  , radii_(config.view<std::vector<double>>(buzzard_sigma_halos::module_label(),
                                            "radii"))
{}

void
buzzard_sigma_halos::set_sample(DataBlock& sample)
{
  // If we had a data member of type std::optional<X>, we would set the
  // value using std::optional::emplace(...) here. emplace takes a set
  // of arguments that it passes to the constructor of X.
  dv_do_dz.emplace(sample);
  hmf.emplace(sample);
  omega_z_sdss.emplace(sample);
  sigma.emplace(sample);
}

std::vector<double>
buzzard_sigma_halos::operator()(double zt, double lnM) const
{
  // For any data members of type std::optional<X>, we have to use operator*
  // to access the X object (as if we were dereferencing a pointer).
  std::vector<double> results(1 + radii_.size());
  double Nval = (*dv_do_dz)(zt) * (*hmf)(lnM, zt) * (*omega_z_sdss)(zt);
  // Number counts followed by the radius bins, repeating over zo bins
  for (std::size_t j = 0; j != radii_.size(); j++) {
    results[j + 1] = (*sigma)(radii_[j], lnM, zt) * Nval;
  }
  results[0] = Nval; // results[i*(radii_.size()+1)+1];
  return results;
}

//
void
buzzard_sigma_halos::finalize_sample(
  cosmosis::DataBlock& sample,
  std::vector<integration_results_v> const& results) const
{
  std::vector<int> NM_status;
  std::vector<int> NM_nregions;
  std::vector<int> NM_nevals;
  std::vector<double> N_vals;
  std::vector<double> N_errors;
  std::vector<double> N_probs;

  std::vector<double> totSigma_vals_temp;
  std::vector<double> totSigma_errors_temp;
  std::vector<double> totSigma_probs_temp;

  std::size_t n_radii_bins = radii_.size();
  for (auto const& result : results) {
    NM_status.push_back(result.status);
    NM_nregions.push_back(result.nregions);
    NM_nevals.push_back(result.neval);

    N_vals.push_back(result.value[0]);
    N_errors.push_back(result.error[0]);
    N_probs.push_back(result.prob[0]);

    totSigma_vals_temp.insert(totSigma_vals_temp.end(),
                              result.value.begin() + 1,
                              result.value.begin() + (n_radii_bins + 1));
    totSigma_errors_temp.insert(totSigma_errors_temp.end(),
                                result.error.begin() + 1,
                                result.error.begin() + (n_radii_bins + 1));
    totSigma_probs_temp.insert(totSigma_probs_temp.end(),
                               result.prob.begin() + 1,
                               result.prob.begin() + (n_radii_bins + 1));
  }
  std::vector<std::size_t> extents{results.size(), n_radii_bins};
  cosmosis::ndarray<double> totSigma_vals(totSigma_vals_temp, extents);
  cosmosis::ndarray<double> totSigma_errors(totSigma_errors_temp, extents);
  cosmosis::ndarray<double> totSigma_probs(totSigma_probs_temp, extents);

  sample.put_val(module_label(), "N_vals", N_vals);
  sample.put_val(module_label(), "N_errors", N_errors);
  sample.put_val(module_label(), "N_probs", N_probs);
  sample.put_val(module_label(), "NM_status", NM_status);
  sample.put_val(module_label(), "NM_nregions", NM_nregions);
  sample.put_val(module_label(), "NM_nevals", NM_nevals);
  sample.put_val(module_label(), "Sigma_radius", radii_);

  sample.put_val(module_label(), "totSigma_vals", totSigma_vals);
  sample.put_val(module_label(), "totSigma_errors", totSigma_errors);
  sample.put_val(module_label(), "totSigma_probs", totSigma_probs);
}

char const*
buzzard_sigma_halos::module_label()
{
  return "buzzard_sigma_halos";
}

// The implementation of make_integration_volumes can be almost the same for
// any CosmoSISIntegrand-type class. Only the names and number of the parameters
// provided need to be changed. It is critical that the names be given in the
// order that correspond to the order of arguments in the class's function call
// operator. While the compiler can verify the number of arguments provided is
// correct, it can not verify that their order matches the order of arguments in
// the function call operator.
std::vector<buzzard_sigma_halos::volume_t>
buzzard_sigma_halos::make_integration_volumes(cosmosis::DataBlock& cfg)
{
  try {
    return y3_cluster::make_integration_volumes_wall_of_numbers(
      cfg, module_label(), "zt", "lnm");
  }
  catch (std::exception const& ex) {
    std::cerr << ex.what();
    throw;
  };
}

DEFINE_COSMOSIS_VECTOR_INTEGRATION_MODULE(buzzard_sigma_halos);
