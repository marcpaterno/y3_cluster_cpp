#ifndef Y3CLUSTER_CLUSTERS_MODULE_HH
#define Y3CLUSTER_CLUSTERS_MODULE_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "/cosmosis/cubacpp/cubacpp/cubacpp.hh"

#include "datablock_reader.hh"
#include "interp_1d.hh"
#include "fits_loader.hh"
#include "gamma_t.hh"

#include <memory>

namespace y3_cluster {
  template <class MODELS>
  class ClustersModule {
    // The radii from the cluster center to predict gamma_t values for
    std::vector<double> radius_bins;
    // The distributions of weak lensing source bins
    std::vector<std::shared_ptr<Interp1D const>> pzsources;
    // The cluster bins in richness-space
    std::vector<y3_cluster::IntegrationRange> lo_bins;
    // The cluster bins in redshift-space
    std::vector<y3_cluster::IntegrationRange> zo_bins;

  public:
    explicit ClustersModule(cosmosis::DataBlock& config);
    void execute(cosmosis::DataBlock& sample);
  };
}


namespace {
  inline std::vector<y3_cluster::IntegrationRange>
  into_ranges(cosmosis::DataBlock db, std::string name)
  {
    auto const edges = get_datablock<std::vector<double>>(db, OPTION_SECTION,
                                                          (name + "_bins").c_str());

    auto ret = std::vector<y3_cluster::IntegrationRange> (edges.size () - 1);

    std::transform  (begin (edges),  end (edges) - 1,
                     begin (edges) + 1,
                     begin (ret),
                     [] (double const &a,  double const &b)
                     {  return y3_cluster::IntegrationRange {a, b}; });

    return ret;
  }

  inline std::vector<double>
  into_edges(const std::vector<y3_cluster::IntegrationRange>& ranges)
  {
    std::vector<double> output;
    for (const auto& ir : ranges)
      output.push_back(ir.transform(0.0));

    output.push_back(ranges[ranges.size() - 1].transform(1.0));
    return output;
  }
}

template <class MODELS>
y3_cluster::ClustersModule<MODELS>::ClustersModule(cosmosis::DataBlock& config)
  // TODO: Possibly set up any optional parameters, like integration params?
  : pzsources(load_pzsources(
              get_datablock<std::string>(config,
                                         "cluster_abundance",
                                         "y3_observables")))
  , lo_bins(into_ranges(config, "lo"))
  , zo_bins(into_ranges(config, "zo"))
{
  auto const radii = get_datablock<std::vector<double>>(config, OPTION_SECTION, "radius_bins");

  radius_bins = std::vector<double>(radii.begin(), radii.end());
}

// TODO:
// - Combine centered and miscentered terms
// - Compute covariance matrix
// - Compute likelihood - done in Michel's module
//
// Then you've got science!
template <class MODELS>
void
y3_cluster::ClustersModule<MODELS>::execute(cosmosis::DataBlock& sample)
{
  // First - output the important parameters
  sample.put_val<std::vector<double>>("cluster_abundance", "lo_bin_edges", into_edges(lo_bins));
  sample.put_val<std::vector<double>>("cluster_abundance", "zo_bin_edges", into_edges(zo_bins));
  sample.put_val<std::vector<double>>("cluster_abundance", "radius_bins", radius_bins);
  sample.put_val<int>("cluster_abundance", "npzsources", pzsources.size());

  // FIXME: Just a test placeholder! Should these come from CosmoSIS?
  double const epsrel = 1.0e-3;
  double const epsabs = 1.0e-12;
  // We have evdidence that at 64 mpi ranks (nut not at 32!) that cuba is forking. so, prepare to do so.
  cubacores(0,0);
  cubacpp::Cuhre c;
  c.maxeval = 100000000;

  // Initialize our big computation stuff up here, so any DataBlock access
  // errors happen up front
  auto integrand = Gamma_T_Integrand<MODELS>
                     {sample, radius_bins, pzsources, lo_bins, zo_bins};

  // Compute abundance counts and gamma_t's
  auto centered_result = integrand.integrate_centered(c, epsrel, epsabs);
  auto miscentered_result = integrand.integrate_miscentered(c, epsrel, epsabs);

  if (centered_result.status != 0)
    throw std::runtime_error("Centered integration did not converge!");
  if (miscentered_result.status != 0)
    throw std::runtime_error("Miscentered integration did not converge!");

  std::cout << "Centered:\n" << centered_result;
  std::cout << "Miscentered:\n" << miscentered_result;

  // All of our output vectors, centered and miscentered for each of:
  //  - gamma_ts
  //  - Cluster abundance predictions (counts)
  //  - Cluster abundance integration error
  //  - Biased cluster abundance predictions (Nb)
  std::vector<double> centered_gamma_ts,
                      centered_cluster_counts,
                      centered_count_errors,
                      centered_biased_counts,
                      miscentered_gamma_ts,
                      miscentered_cluster_counts,
                      miscentered_count_errors,
                      miscentered_biased_counts;

  // Sort centered abundance counts/biased counts/error bars and gamma_t's
  for (const auto& bin : centered_result) {
    centered_cluster_counts.push_back(bin.N);
    centered_count_errors.push_back(bin.N_error);
    centered_biased_counts.push_back(bin.Nb);
    for (const auto& zsrc_bin : bin.gamma_ts)
      for (const auto gt : zsrc_bin)
        centered_gamma_ts.push_back(gt);
  }

  // Sort miscentered abundance counts/biased counts/error bars and gamma_t's
  for (const auto& bin : miscentered_result) {
    miscentered_cluster_counts.push_back(bin.N);
    miscentered_count_errors.push_back(bin.N_error);
    miscentered_biased_counts.push_back(bin.Nb);
    for (const auto& zsrc_bin : bin.gamma_ts)
      for (const auto gt : zsrc_bin)
        miscentered_gamma_ts.push_back(gt);
  }

  // Store abundance counts and gamma_t's
  sample.put_val<std::vector<double>>("cluster_abundance", "centered_gamma_ts",
                                      centered_gamma_ts);
  sample.put_val<std::vector<double>>("cluster_abundance", "centered_cluster_counts",
                                      centered_cluster_counts);
  sample.put_val<std::vector<double>>("cluster_abundance", "centered_cluster_count_errors",
                                      centered_count_errors);
  sample.put_val<std::vector<double>>("cluster_abundance", "centered_cluster_biased_counts",
                                      centered_biased_counts);
  sample.put_val<std::vector<double>>("cluster_abundance", "miscentered_gamma_ts",
                                      miscentered_gamma_ts);
  sample.put_val<std::vector<double>>("cluster_abundance", "miscentered_cluster_counts",
                                      miscentered_cluster_counts);
  sample.put_val<std::vector<double>>("cluster_abundance", "miscentered_cluster_count_errors",
                                      miscentered_count_errors);
  sample.put_val<std::vector<double>>("cluster_abundance", "miscentered_cluster_biased_counts",
                                      miscentered_biased_counts);
}

#endif
