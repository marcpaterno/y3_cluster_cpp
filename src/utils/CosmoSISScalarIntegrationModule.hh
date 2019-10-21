#ifndef Y3_CLUSTER_COSMOSIS_SCALAR_INTEGRATION_MODULE_HH
#define Y3_CLUSTER_COSMOSIS_SCALAR_INTEGRATION_MODULE_HH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/entry.hh"
#include "cubacpp/cuhre.hh"
#include "cubacpp/vegas.hh"
#include "cubacpp/suave.hh"
#include "cubacpp/integrand_traits.hh"
#include "cubacpp/integration_result.hh"
#include "cubacpp/integration_volume.hh"

#include <iostream>
#include <vector>

namespace y3_cluster {
  // CosmoSISScalarIntegrationModule is a class template used to create a class
  // that is a CosmoSIS physics module, and which, for each CosmoSIS sample,
  // calculates the definite integral of a scalar-valued function.
  //
  // A module instantiated from this class template can be configured to
  // integrate multiple distinct volumes (ranges of integration for each
  // variable of integration). In addition, it can be configured to do these
  // integrations on each point in a multidimensional grid.
  //
  // The results are put into the cosmosis::DataBlock for the sample
  //
  // Likelihood calculations based on the result of this integration must be
  // done by a separate CosmoSIS module.

  template <typename COSMOSISINTEGRAND, typename ALG = cubacpp::Cuhre>
  class CosmoSISScalarIntegrationModule {
  public:
    using IntegrandType = COSMOSISINTEGRAND;
    using Algorithm = ALG;
    using volume_t = cubacpp::integration_volume_for_t<IntegrandType>;
    using grid_point_t = typename IntegrandType::grid_point_t;
    using integration_results_t = cubacpp::integration_result;

    explicit CosmoSISScalarIntegrationModule(cosmosis::DataBlock& cfg);

    void execute(cosmosis::DataBlock& sample);

  private:
    void integrate_one_volume(volume_t const& vol,
                              std::vector<integration_results_t>& results);

    // finalize_sample() is where products can be put into the
    // cosmosis::DataBlock representing the current sample. The object 'sample'
    // passed to this function will be the exact same object as was passed to
    // the most recent call to set_sample(). The object 'results' will be the
    // results of the integration that has just been done for that sample. This
    // is generally the item which should be put into the sample.
    void finalize_sample(
      cosmosis::DataBlock& sample,
      std::vector<cubacpp::integration_result> const& results) const;

    IntegrandType integrand_;
    Algorithm algorithm_;
    std::vector<volume_t> volumes_;
    std::vector<grid_point_t> grid_points_;
    double eps_rel_;
    double eps_abs_;
  };
} // namespace y3_cluster

template <typename I, typename A>
y3_cluster::CosmoSISScalarIntegrationModule<I, A>::
  CosmoSISScalarIntegrationModule(cosmosis::DataBlock& cfg)
try : integrand_(cfg),
      algorithm_(),
      volumes_(IntegrandType::make_integration_volumes(cfg)),
      grid_points_(IntegrandType::make_grid_points(cfg)),
      eps_rel_(cfg.view<double>(IntegrandType::module_label(), "eps_rel")),
      eps_abs_(cfg.view<double>(IntegrandType::module_label(), "eps_abs")) {
  algorithm_.maxeval = cfg.view<int>(IntegrandType::module_label(), "max_eval");
  cubacores(0, 0);
}
catch (cosmosis::Exception const&) {
  std::cerr
    << "\nDuring construction of a CosmoSISScalarIntegrationModule, the "
       "lookup of some parameter"
    << "\nfailed. It may be a wrong name, or a wrong type.\n";
}
catch (...) {
  std::cerr << "\nUnknown exception type thrown while constructing a "
               "CosmoSISScalarIntegrationModule.\n\n";
}

template <typename I, typename A>
void
y3_cluster::CosmoSISScalarIntegrationModule<I, A>::integrate_one_volume(
  volume_t const& volume,
  std::vector<integration_results_t>& results)
{
  for (auto const& grid_point : grid_points_) {
    integrand_.set_grid_point(grid_point);
    auto result = algorithm_.integrate(integrand_, eps_rel_, eps_abs_, volume);
    results.push_back(result);
  }
}

template <typename I, typename A>
void
y3_cluster::CosmoSISScalarIntegrationModule<I, A>::execute(
  cosmosis::DataBlock& sample)
{
  // Prepare the integrand for this sample.
  integrand_.set_sample(sample);

  std::vector<integration_results_t> results;
  results.reserve(volumes_.size() * grid_points_.size());

  for (auto const& volume : volumes_) {
    integrate_one_volume(volume, results);
  }

  // Put the results into the sample.
  finalize_sample(sample, results);
}

template <typename I, typename A>
void
y3_cluster::CosmoSISScalarIntegrationModule<I, A>::finalize_sample(
  cosmosis::DataBlock& sample,
  std::vector<cubacpp::integration_result> const& res) const
{
  using cosmosis::ndarray;
  using cubacpp::integration_result;

  auto ngrid_points = grid_points_.size();
  auto nvolumes = volumes_.size();
  auto nresults = nvolumes * ngrid_points;

  // Create ndarray to give a view into the 'res' vector.
  ndarray<integration_result> results(res, {nvolumes, ngrid_points});

  // Create storage buffers for ndarrays to be inserted into the sample, and
  // then create the ndarrays.
  std::vector<double> vals_buffer(nresults);
  ndarray<double> vals(vals_buffer, {nvolumes, ngrid_points});

  std::vector<double> errors_buffer(nresults);
  ndarray<double> errors(errors_buffer, {nvolumes, ngrid_points});

  std::vector<double> probs_buffer(nresults);
  ndarray<double> probs(probs_buffer, {nvolumes, ngrid_points});

  std::vector<int> statuses_buffer(nresults);
  ndarray<int> statuses(statuses_buffer, {nvolumes, ngrid_points});

  std::vector<int> nregions_buffer(nresults);
  ndarray<int> nregions(nregions_buffer, {nvolumes, ngrid_points});

  for (std::size_t ivol = 0; ivol != nvolumes; ++ivol) {
    for (std::size_t jgp = 0; jgp != ngrid_points; ++jgp) {
      vals(ivol, jgp) = results(ivol, jgp).value;
      errors(ivol, jgp) = results(ivol, jgp).error;
      probs(ivol, jgp) = results(ivol, jgp).prob;
      statuses(ivol, jgp) = results(ivol, jgp).status;
      nregions(ivol, jgp) = results(ivol, jgp).nregions;
    }
  }

  auto module_label = IntegrandType::module_label();
  sample.put_val(module_label, "vals", vals);
  sample.put_val(module_label, "errors", errors);
  sample.put_val(module_label, "probs", probs);
  sample.put_val(module_label, "status", statuses);
  sample.put_val(module_label, "nregions", nregions);
}

#endif
