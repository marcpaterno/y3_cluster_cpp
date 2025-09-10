#include "utils/make_integration_volumes.hh"
#include "utils/module_macros.hh"

#include "cosmosis/datablock/datablock.hh"
#include "cubacpp/integration_result.hh"
#include "cubacpp/integration_volume.hh"

#include "models/dv_do_dz_t.hh"
#include "models/hmf_t.hh"
#include "models/int_lc_lt_des_t.hh"
#include "models/int_zo_zt_des_t.hh"
#include "models/mor_des_2022_logm.hh"
#include "models/omega_z_des.hh"
#include "utils/make_grid_points.hh"
#include <optional>
#include <vector>

// y1_analysis_logm_2022 is a class that models the concept of "CosmoSISIntegrand",
// and is thus suitable for use as the template parameter for the class template
// CosmosisIntegrationModule.
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
class y1_analysis_logm_2022 {
private:
  // We define the type alias volume_t to be the right dimensionality
  // of integration volume for our integrand. If we were to change the
  // number of arguments required by the function call operator (below),
  // we would need to also modify this type alias to keep consistent.
  using volume_t = cubacpp::IntegrationVolume<4>;

  // State obtained from each sample.
  // If there were a type X that did not have a default constructor,
  // we would use std::optional<X> as our data member.
  std::optional<y3_cluster::INT_LC_LT_DES_t> lc_lt;
  std::optional<y3_cluster::MOR_DES_2022_LOGM> mor;
  std::optional<y3_cluster::OMEGA_Z_DES> omega_z;
  std::optional<y3_cluster::DV_DO_DZ_t> dv_do_dz;
  std::optional<y3_cluster::HMF_t> hmf;
  std::optional<y3_cluster::INT_ZO_ZT_DES_t> int_zo_zt;
  std::vector<double> zo_low_;
  std::vector<double> zo_high_;

public:
  using grid_t = y3_cluster::grid_t<1>;
  using grid_point_t = grid_t::value_type;

  // Initialize my integrand object from the parameters read
  // from the relevant block in the CosmoSIS ini file.
  explicit y1_analysis_logm_2022(cosmosis::DataBlock& config);

  // Set any data members from values read from the current sample.
  // Do not attempt to copy the sample!.
  void set_sample(cosmosis::DataBlock& sample);

  // The function to be integrated. All arguments to this function must be of
  // type double, and there must be at least two of them (because our
  // integration routine does not work for functions of one variable). The
  // function is const because calling it does not change the state of the
  // object.
  std::vector<double> operator()(double lo,
                                 double lt,
                                 double zt,
                                 double lnM) const;

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

y1_analysis_logm_2022::y1_analysis_logm_2022(DataBlock& config)
  : lc_lt()
  , mor()
  , omega_z()
  , dv_do_dz()
  , hmf()
  , int_zo_zt()
  , zo_low_(
      config.view<std::vector<double>>(y1_analysis_logm_2022::module_label(), "zo_low"))
  , zo_high_(
      config.view<std::vector<double>>(y1_analysis_logm_2022::module_label(), "zo_high"))
{}

void
y1_analysis_logm_2022::set_sample(DataBlock& sample)
{
  // If we had a data member of type std::optional<X>, we would set the
  // value using std::optional::emplace(...) here. emplace takes a set
  // of arguments that it passes to the constructor of X.
  lc_lt.emplace(sample);
  mor.emplace(sample);
  dv_do_dz.emplace(sample);
  hmf.emplace(sample);
  omega_z.emplace(sample);
}

std::vector<double>
y1_analysis_logm_2022::operator()(double lo, double lt, double zt, double lnM) const
{
  // For any data members of type std::optional<X>, we have to use operator*
  // to access the X object (as if we were dereferencing a pointer).
  std::vector<double> results(2 * zo_low_.size());
  double mass = std::exp(lnM);
  double common_term = (*lc_lt)(lo, lt, zt) * (*mor)(lt, lnM, zt) *
                       (*dv_do_dz)(zt) * (*hmf)(lnM, zt) * (*omega_z)(zt);
  // Number counts first in the returned results, followed by the masses
  for (std::size_t i = 0; i != zo_low_.size(); i++) {
    results[i] = (*int_zo_zt)(zo_low_[i], zo_high_[i], zt) * common_term;
    results[i + zo_low_.size()] = mass * results[i];
  }
  return results;
}

//
void
y1_analysis_logm_2022::finalize_sample(
  cosmosis::DataBlock& sample,
  std::vector<integration_results_v> const& results) const
{
  std::vector<int> NM_status;
  std::vector<int> NM_nregions;
  std::vector<int> NM_nevals;

  std::vector<double> N_vals_temp;
  std::vector<double> N_errors_temp;
  std::vector<double> N_probs_temp;
  std::vector<double> totM_vals_temp;
  std::vector<double> totM_errors_temp;
  std::vector<double> totM_probs_temp;

  ///// how do I know how many dimensions zo_low_ or zo_high_ or radiii have?
  //
  // TODO: Try to refactor this code, to abstract away the manipulations done
  // for each physics quantity.
  std::size_t n_zo_bins = zo_low_.size();
  for (auto const& result : results) {
    NM_status.push_back(result.status);
    NM_nregions.push_back(result.nregions);
    NM_nevals.push_back(result.neval);

    N_vals_temp.insert(N_vals_temp.end(),
                       result.value.begin(),
                       result.value.begin() + n_zo_bins);
    N_errors_temp.insert(N_errors_temp.end(),
                         result.error.begin(),
                         result.error.begin() + n_zo_bins);
    N_probs_temp.insert(
      N_probs_temp.end(), result.prob.begin(), result.prob.begin() + n_zo_bins);

    totM_vals_temp.insert(totM_vals_temp.end(),
                          result.value.begin() + n_zo_bins,
                          result.value.end());
    totM_errors_temp.insert(totM_errors_temp.end(),
                            result.error.begin() + n_zo_bins,
                            result.error.end());
    totM_probs_temp.insert(totM_probs_temp.end(),
                           result.prob.begin() + n_zo_bins,
                           result.prob.end());
  }
  std::vector<std::size_t> extents{results.size(), n_zo_bins};
  cosmosis::ndarray<double> N_vals(N_vals_temp, extents);
  cosmosis::ndarray<double> N_errors(N_errors_temp, extents);
  cosmosis::ndarray<double> N_probs(N_probs_temp, extents);
  cosmosis::ndarray<double> totM_vals(totM_vals_temp, extents);
  cosmosis::ndarray<double> totM_errors(totM_errors_temp, extents);
  cosmosis::ndarray<double> totM_probs(totM_probs_temp, extents);

  sample.put_val(module_label(), "N_vals", N_vals);
  sample.put_val(module_label(), "N_errors", N_errors);
  sample.put_val(module_label(), "N_probs", N_probs);
  sample.put_val(module_label(), "NM_status", NM_status);
  sample.put_val(module_label(), "NM_nregions", NM_nregions);
  sample.put_val(module_label(), "NM_nevals", NM_nevals);

  sample.put_val(module_label(), "totM_vals", totM_vals);
  sample.put_val(module_label(), "totM_errors", totM_errors);
  sample.put_val(module_label(), "totM_probs", totM_probs);
}

char const*
y1_analysis_logm_2022::module_label()
{
  return "y1_analysis_logm_2022";
}

// The implementation of make_integration_volumes can be almost the same for
// any CosmoSISIntegrand-type class. Only the names and number of the parameters
// provided need to be changed. It is critical that the names be given in the
// order that correspond to the order of arguments in the class's function call
// operator. While the compiler can verify the number of arguments provided is
// correct, it can not verify that their order matches the order of arguments in
// the function call operator.
std::vector<y1_analysis_logm_2022::volume_t>
y1_analysis_logm_2022::make_integration_volumes(cosmosis::DataBlock& cfg)
{
  try {
    return y3_cluster::make_integration_volumes_wall_of_numbers(
      cfg, module_label(), "lo", "lt", "zt", "lnm");
  }
  catch (std::exception const& ex) {
    std::cerr << ex.what();
    throw;
  };
}

DEFINE_COSMOSIS_VECTOR_INTEGRATION_MODULE(y1_analysis_logm_2022)
