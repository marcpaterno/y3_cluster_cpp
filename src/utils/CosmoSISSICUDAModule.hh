#ifndef Y3_CLUSTER_COSMOSIS_SI_CUDA_MODULE_HH
#define Y3_CLUSTER_COSMOSIS_SI_CUDA_MODULE_HH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/entry.hh"
#include "cosmosis/datablock/ndarray.hh"

#include "cubacpp/arity.hh"
#include "cubacpp/integration_result.hh"

#include "utils/make_grid_points.hh"
#include "utils/make_cuda_integration_volumes.cuh"
#include "utils/gpu_integrator.cuh"

#include "cudaPagani/quad/util/Volume.cuh"
#include "cudaPagani/quad/GPUquad/Pagani.cuh"
#include "vegas/vegasT.cuh"

#include <iostream>
#include <stdexcept>
#include <vector>


namespace y3_cluster {
  // CosmoSISSICUDAModule is a class template used to create a class
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

  template <typename IGRAND>
  class CosmoSISSICUDAModule {
  public:
    using IntegrandType = IGRAND;
    static constexpr int ndim = cubacpp::arity<IntegrandType>();
    using volume_t = quad::Volume<double, ndim>;
    using grid_point_t = typename IntegrandType::grid_point_t;
    using my_grid_t = grid_t<std::tuple_size<grid_point_t>::value>;
    using integration_results_t = cubacpp::integration_result;

    explicit CosmoSISSICUDAModule(cosmosis::DataBlock& cfg);

    // Evaluate all the integrals specified (all grid points, and all volumes),
    // and record the results in the sample.
    void execute(cosmosis::DataBlock& sample);

  private:
    // Return the number of results that will be entered into the DataBlock.
    // This is the total number of integrals computed in one call to execute().
    std::size_t num_results() const;

    std::vector<integration_results_t>
    integrate_cartesian_product_of_volumes_and_gridpoints();

    std::vector<integration_results_t>
    integrate_zipped_sequence_of_volumes_and_gridpoints();

    // finalize_sample() is where products can be put into the
    // cosmosis::DataBlock representing the current sample. The object 'sample'
    // passed to this function will be the exact same object as was passed to
    // the most recent call to set_sample(). The object 'results' will be the
    // results of the integration that has just been done for that sample. This
    // is generally the item which should be put into the sample.
    void finalize_sample(
      cosmosis::DataBlock& sample,
      std::vector<integration_results_t> const& results) const;

    void finalize_sample_cartesian_product_of_volumes_and_gridpoints(
      cosmosis::DataBlock& sample,
      std::vector<integration_results_t> const& results,
      my_grid_t const* pgrid) const;

    void
    finalize_sample_cartesian_product_of_volumes_and_gridpoints(
      cosmosis::DataBlock& sample,
      std::vector<integration_results_t> const& results) const
    {
      finalize_sample_cartesian_product_of_volumes_and_gridpoints(
        sample, results, &grid_points_);
    }

    void finalize_sample_zipped_sequence_of_volumes_and_gridpoints(
      cosmosis::DataBlock& sample,
      std::vector<integration_results_t> const& results) const;

    IntegrandType integrand_;
    y3_cuda::MultiDimensionalIntegrator<ndim> algorithm_;
    std::vector<volume_t> volumes_;
    my_grid_t grid_points_;
    double eps_rel_;
    double eps_abs_;
    bool use_cartesian_product_of_volumes_and_gridpoints_;
  };
} // namespace y3_cluster

template <typename I>
y3_cluster::CosmoSISSICUDAModule<I>::CosmoSISSICUDAModule(
  cosmosis::DataBlock& cfg)
try : integrand_(cfg),
    // The c'tor for algorithm_ will get adjusted at some later date...
	algorithm_(cfg.view<std::string>(IntegrandType::module_label(), "algorithm")),
  volumes_(IntegrandType::make_integration_volumes(cfg)),
  grid_points_(IntegrandType::make_grid_points(cfg)),
  eps_rel_(cfg.view<double>(IntegrandType::module_label(), "eps_rel")),
  eps_abs_(cfg.view<double>(IntegrandType::module_label(), "eps_abs")),
  use_cartesian_product_of_volumes_and_gridpoints_(
    cfg.view<bool>(IntegrandType::module_label(), "use_cartesian_product")) {
  if (not use_cartesian_product_of_volumes_and_gridpoints_) {
    if (volumes_.size() != grid_points_.size()) {
      throw std::runtime_error(
        "An integration module was configured to use a zipped sequence of "
        "volumes and gridpoints,\n"
        "but the number of volumes did not equal the number of gridpoints.\n");
    }
  }
  algorithm_.set_maxeval(cfg.view<double>(IntegrandType::module_label(), "max_eval"));
}
catch (cosmosis::Exception const&) {
  std::cerr
    << "\nDuring construction of a CosmoSISSICUDAModule, the "
       "lookup of some parameter"
    << "\nfailed. It may be a wrong name, or a wrong type.\n";
}
catch (std::exception const& e) {
  std::cerr << "\nDuring construction of a CosmoSISSICUDAModule, an "
               "std::exeption throw was encountered.\n"
               "The error message was:\n"
            << e.what();
}
catch (...) {
  std::cerr << "\nUnknown exception type thrown while constructing a "
               "CosmoSISSICUDAModule.\n\n";
}

template <typename I>
void
y3_cluster::CosmoSISSICUDAModule<I>::execute(
  cosmosis::DataBlock& sample)
{
  integrand_.set_sample(sample);
  auto results = use_cartesian_product_of_volumes_and_gridpoints_ ?
                   integrate_cartesian_product_of_volumes_and_gridpoints() :
                   integrate_zipped_sequence_of_volumes_and_gridpoints();
  finalize_sample(sample, results);
}

template <typename I>
std::size_t
y3_cluster::CosmoSISSICUDAModule<I>::num_results() const
{
  return (use_cartesian_product_of_volumes_and_gridpoints_) ?
           volumes_.size() * grid_points_.size() :
           volumes_.size();
}

template <typename I>
std::vector<cubacpp::integration_result>
y3_cluster::CosmoSISSICUDAModule<
  I>::integrate_cartesian_product_of_volumes_and_gridpoints()
{
  std::vector<integration_results_t> results;
  results.reserve(num_results());

  for (auto& volume : volumes_) {
    for (auto const& grid_point : grid_points_.points) {
      integrand_.set_grid_point(grid_point);
      cuhreResult res = 
        algorithm_.integrate(integrand_, eps_abs_, eps_rel_, &volume);
      results.push_back(cubacpp::integration_result(
            res.estimate,
            res.errorest,
            0.0,    // probability not returned by PAGANI
            static_cast<long long>(res.neval),
            static_cast<int>(res.nregions),
            static_cast<int>(res.status)));
    }
  }
  return results;
}

template <typename IGRAND>
std::vector<cubacpp::integration_result>
y3_cluster::CosmoSISSICUDAModule<
  IGRAND>::integrate_zipped_sequence_of_volumes_and_gridpoints()
{
  std::vector<integration_results_t> results;
  results.reserve(num_results());
  for (std::size_t i = 0; i != num_results(); ++i) {
    integrand_.set_grid_point(grid_points_.points[i]);
    cuhreResult res = 
      algorithm_.integrate(integrand_, eps_abs_, eps_rel_, &volumes_[i]);
     results.push_back(cubacpp::integration_result(
           res.estimate,
           res.errorest,
           0.0, // probability not returned by PAGANI
           static_cast<long long>(res.neval),
           static_cast<int>(res.nregions),
           static_cast<int>(res.status)));
  }
  return results;
}

template <typename I>
void
y3_cluster::CosmoSISSICUDAModule<I>::finalize_sample(
  cosmosis::DataBlock& sample,
  std::vector<cubacpp::integration_result> const& res) const
{
  std::cerr << "Finalizing sample.\n";
  if (use_cartesian_product_of_volumes_and_gridpoints_)
    finalize_sample_cartesian_product_of_volumes_and_gridpoints(sample, res);
  else
    finalize_sample_zipped_sequence_of_volumes_and_gridpoints(sample, res);
}

template <typename I>
void
y3_cluster::CosmoSISSICUDAModule<I>::
  finalize_sample_zipped_sequence_of_volumes_and_gridpoints(
    cosmosis::DataBlock& sample,
    std::vector<cubacpp::integration_result> const& res) const
{
  using cosmosis::ndarray;
  using cubacpp::integration_result;

  auto const nresults = num_results();

  // Create an ndarray to give a view into the 'res' vector.
  ndarray<integration_result> results(res, {nresults});

  // Create storage buffers for ndarrays to be inserted into the sample, and
  // then create the ndarrays.
  std::vector<double> vals_buffer(nresults);
  ndarray<double> vals(vals_buffer, {nresults});

  std::vector<double> errors_buffer(nresults);
  ndarray<double> errors(errors_buffer, {nresults});

  std::vector<double> probs_buffer(nresults);
  ndarray<double> probs(probs_buffer, {nresults});

  std::vector<int> statuses_buffer(nresults);
  ndarray<int> statuses(statuses_buffer, {nresults});

  std::vector<int> nregions_buffer(nresults);
  ndarray<int> nregions(nregions_buffer, {nresults});

  for (std::size_t i = 0; i != nresults; ++i) {
    vals(i) = results(i).value;
    errors(i) = results(i).error;
    probs(i) = results(i).prob;
    statuses(i) = results(i).status;
    nregions(i) = results(i).nregions;
  }

  auto module_label = IntegrandType::module_label();
  sample.put_val(module_label, "vals", vals);
  sample.put_val(module_label, "errors", errors);
  sample.put_val(module_label, "probs", probs);
  sample.put_val(module_label, "status", statuses);
  sample.put_val(module_label, "nregions", nregions);
}
template <typename I>
void
y3_cluster::CosmoSISSICUDAModule<I>::
  finalize_sample_cartesian_product_of_volumes_and_gridpoints(
    cosmosis::DataBlock& sample,
    std::vector<cubacpp::integration_result> const& res,
    my_grid_t const* pgrid) const
{
  using cosmosis::ndarray;
  using cubacpp::integration_result;

  auto const nresults = num_results();

  auto const ngrid_points = grid_points_.size();
  auto const nvolumes = volumes_.size();

  // Create an ndarray to give a view into the 'res' vector.
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

  std::vector<int> vols_buffer(nresults);
  ndarray<int> vols(vols_buffer, {nvolumes, ngrid_points});

  constexpr std::size_t naxes = std::tuple_size<grid_point_t>::value;
  std::vector<double> gridpoints_buffer(ngrid_points * naxes);
  ndarray<double> gridpoints(gridpoints_buffer, {ngrid_points, naxes});

  if (pgrid) {
    for (std::size_t jgp = 0; jgp != ngrid_points; ++jgp) {
      for (std::size_t iaxis = 0; iaxis != naxes; ++iaxis) {
        gridpoints(jgp, iaxis) = pgrid->points[jgp][iaxis];
      }
    }
  }

  for (std::size_t ivol = 0; ivol != nvolumes; ++ivol) {
    for (std::size_t jgp = 0; jgp != ngrid_points; ++jgp) {
      vols(ivol, jgp) = ivol;
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
  sample.put_val(module_label, "ivol", vols);
  sample.put_val(module_label, "gridpoints", gridpoints);
}

#endif
