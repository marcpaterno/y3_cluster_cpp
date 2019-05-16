#include "module_macros.hh"
#include "make_integration_volumes.hh"

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "cubacpp/integration_result.hh"
#include "cubacpp/integration_volume.hh"

#include "lc_lt_t.hh"
#include "mor_t2.hh"
#include "hmf_t.hh"
#include "dv_do_dz_t.hh"
#include "omega_z_sdss.hh"
#include "int_zo_zt_t.hh"
#include <optional>
#include <vector>
using namespace y3_cluster;

// y1_analysis is a class that models the concept of "CosmoSISIntegrand",
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
class y1_analysis {
private:
  // We define the type alias volume_t to be the right dimensionality
  // of integration volume for our integrand. If we were to change the
  // number of arguments required by the function call operator (below),
  // we would need to also modify this type alias to keep consistent.
  using volume_t = cubacpp::IntegrationVolume<4>;

  // State obtained from each sample.
  // If there were a type X that did not have a default constructor,
  // we would use std::optional<X> as our data member.
  std::optional<LC_LT_t> lc_lt;
  std::optional<MOR_t2> mor;
  std::optional<OMEGA_Z_SDSS> omega_z_sdss;
  std::optional<DV_DO_DZ_t> dv_do_dz;
  std::optional<HMF_t> hmf;
  std::optional<INT_ZO_ZT_t> int_zo_zt;
  std::vector<double> zo_low_;
  std::vector<double> zo_high_;

public:
  // Initialize my integrand object from the parameters read
  // from the relevant block in the CosmoSIS ini file.
  explicit y1_analysis(cosmosis::DataBlock& config);

  // Set any data members from values read from the current sample.
  // Do not attempt to copy the sample!.
  void set_sample(cosmosis::DataBlock& sample);

  // The function to be integrated. All arguments to this function must be of
  // type double, and there must be at least two of them (because our
  // integration routine does not work for functions of one variable). The
  // function is const because calling it does not change the state of the
  // object.
  std::vector<double> operator()(double lo, double lt, double zt, double lnM) const;

  // finalize_sample() is where products can be put into the cosmosis::DataBlock
  // representing the current sample. The object 'sample' passed to this
  // function will be the exact same object as was passed to the most recent
  // call to set_sample(). The object 'results' will be the results of the
  // integration that has just been done for that sample. This is generally the
  // item which should be put into the sample.
  void finalize_sample(
    cosmosis::DataBlock& sample,
    std::vector<cubacpp::integration_results_v> const& results) const;

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

// We put the module label "modulelabel" in an anonymous namespace to make sure
// no other compilation unit can see it, and so that it has static lifetime.
// This means it is created when the library is loaded, and is destroyed when
// the program is shut down.
namespace {
  // We make the modulelabel const because it should never be modified.
  std::string const modulelabel("y1_analysis");
} // anonymous namespace

y1_analysis::y1_analysis(DataBlock& config) : 
	lc_lt(),
       	mor(),
       	omega_z_sdss(),
       	dv_do_dz(),
       	hmf(), 
        int_zo_zt(),
	zo_low_(config.view<std::vector<double>>(modulelabel, "zo_low")),
        zo_high_(config.view<std::vector<double>>(modulelabel, "zo_high")){}

void
y1_analysis::set_sample(DataBlock& sample)
{
  // If we had a data member of type std::optional<X>, we would set the
  // value using std::optional::emplace(...) here. emplace takes a set
  // of arguments that it passes to the constructor of X.
  lc_lt.emplace(sample);
  mor.emplace(sample);
  dv_do_dz.emplace(sample);
  hmf.emplace(sample);
  omega_z_sdss.emplace(sample);
  int_zo_zt.emplace(sample);
}

std::vector<double>
y1_analysis::operator()(double lo, double lt, double zt, double lnM) const
{
  // For any data members of type std::optional<X>, we have to use operator*
  // to access the X object (as if we were dereferencing a pointer).
  std::vector<double> results(2 * zo_low_.size());
  double mass = std::exp(lnM);
  double common_term = (*lc_lt)(lo, lt, zt)
	  	     * (*mor)(lt, lnM, zt)
		     * (*dv_do_dz)(zt)
		     * (*hmf)(lnM, zt)
		     * (*omega_z_sdss)(zt);
  // Number counts first in the returned results, followed by the masses 
  for (std::size_t i =0; i != zo_low_.size(); i++){
  	results[i]= (*int_zo_zt)(zo_low_[i], zo_high_[i], zt)
			* common_term;
        results[i+zo_low_.size()] = mass * results[i];
  }
  return results;
}

//
void
y1_analysis::finalize_sample(
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
  std::size_t n_zo_bins = zo_low_.size();
  for (auto const& result : results) {
    NM_status.push_back(result.status);
    NM_nregions.push_back(result.nregions);
    NM_nevals.push_back(result.neval);

    N_vals_temp.insert(N_vals_temp.end(), result.value.begin(), result.value.begin()+n_zo_bins);
    N_errors_temp.insert(N_errors_temp.end(), result.error.begin(), result.error.begin()+n_zo_bins);
    N_probs_temp.insert(N_probs_temp.end(), result.prob.begin(), result.prob.begin()+n_zo_bins);

    totM_vals_temp.insert(totM_vals_temp.end(), result.value.begin()+n_zo_bins, result.value.end());
    totM_errors_temp.insert(totM_errors_temp.end(), result.error.begin()+n_zo_bins, result.error.end());
    totM_probs_temp.insert(totM_probs_temp.end(), result.prob.begin()+n_zo_bins, result.prob.end());
  }
  std::vector<std::size_t> extents{results.size(), n_zo_bins};
  cosmosis::ndarray<double> N_vals(N_vals_temp, extents);
  cosmosis::ndarray<double> N_errors(N_errors_temp, extents);
  cosmosis::ndarray<double> N_probs(N_probs_temp, extents);
  cosmosis::ndarray<double> totM_vals(totM_vals_temp, extents);
  cosmosis::ndarray<double> totM_errors(totM_errors_temp, extents);
  cosmosis::ndarray<double> totM_probs(totM_probs_temp, extents);

  sample.put_val(modulelabel, "N_vals", N_vals);
  sample.put_val(modulelabel, "N_errors", N_errors);
  sample.put_val(modulelabel, "N_probs", N_probs);
  sample.put_val(modulelabel, "NM_status", NM_status);
  sample.put_val(modulelabel, "NM_nregions", NM_nregions);
  sample.put_val(modulelabel, "NM_nevals", NM_nevals);

  sample.put_val(modulelabel, "totM_vals", totM_vals);
  sample.put_val(modulelabel, "totM_errors", totM_errors);
  sample.put_val(modulelabel, "totM_probs", totM_probs);
}

// The implementation of make_integration_volumes can be almost the same for
// any CosmoSISIntegrand-type class. Only the names and number of the parameters
// provided need to be changed. It is critical that the names be given in the
// order that correspond to the order of arguments in the class's function call
// operator. While the compiler can verify the number of arguments provided is
// correct, it can not verify that their order matches the order of arguments in
// the function call operator.
std::vector<y1_analysis::volume_t>
y1_analysis::make_integration_volumes(cosmosis::DataBlock& cfg)
{
  try {
    return y3_cluster::make_integration_volumes(
      cfg, modulelabel, "lo", "lt", "zt", "lnm");
  }
  catch (std::exception const& ex) {
    std::cerr << ex.what();
    throw;
  };
}

DEFINE_COSMOSIS_MODULE(y1_analysis);
