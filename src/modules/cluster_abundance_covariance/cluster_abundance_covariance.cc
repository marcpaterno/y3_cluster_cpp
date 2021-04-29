#include "cosmosis/datablock/datablock.hh"

#include "models/default_models.hh"
#include "models/sample_variance.hh"
#include "utils/datablock_reader.hh"
#include "utils/integration_range.hh"

using namespace y3_cluster;

namespace y3_cluster {
  class ClusterCovarianceModule {
    std::vector<IntegrationRange> lo_bins;
    std::vector<IntegrationRange> zo_bins;

  public:
    explicit ClusterCovarianceModule(cosmosis::DataBlock& config);
    template <typename OMEGA_Z>
    DATABLOCK_STATUS execute(cosmosis::DataBlock& sample);
  };
}

static inline std::vector<IntegrationRange>
_get_ranges(cosmosis::DataBlock db, std::string name)
{
  auto const edges = get_datablock<std::vector<double>>(
    db, OPTION_SECTION, (name + "_bins").c_str());

  auto ret = std::vector<IntegrationRange>(edges.size() - 1);

  std::transform(begin(edges),
                 end(edges) - 1,
                 begin(edges) + 1,
                 begin(ret),
                 [](double const& a, double const& b) {
                   return y3_cluster::IntegrationRange{a, b};
                 });

  return ret;
}

y3_cluster::ClusterCovarianceModule::ClusterCovarianceModule(
  cosmosis::DataBlock& config)
  : lo_bins(_get_ranges(config, "lo")), zo_bins(_get_ranges(config, "zo"))
{}

template <typename OMEGA_Z>
DATABLOCK_STATUS
y3_cluster::ClusterCovarianceModule::execute(cosmosis::DataBlock& sample)
{
  std::size_t const NRICHNESS = lo_bins.size(), NREDSHIFT = zo_bins.size();

  OMEGA_Z omega_z(sample);
  SampleVariance_t sv(sample, omega_z, zo_bins);

  // Create cluster count covariance - initialized to zeroes
  std::size_t ncounts = NRICHNESS * NREDSHIFT;
  std::vector<double> centered_cluster_covariance(ncounts * ncounts),
    miscentered_cluster_covariance(ncounts * ncounts);
  for (auto i = 0u; i < ncounts * ncounts; i++) {
    centered_cluster_covariance[i] = 0.0;
    miscentered_cluster_covariance[i] = 0.0;
  }

  // Pull output we need:
  //  * (mis)centered cluster abundance prediction
  //  * (mis)centered cluster abundance integration error
  //  * (mis)centered biased cluster abundance prediction
  const auto centered_counts = get_datablock<std::vector<double>>(
               sample, "cluster_abundance", "centered_cluster_counts"),
             miscentered_counts = get_datablock<std::vector<double>>(
               sample, "cluster_abundance", "miscentered_cluster_counts"),
             centered_biased_counts = get_datablock<std::vector<double>>(
               sample, "cluster_abundance", "centered_cluster_biased_counts"),
             miscentered_biased_counts = get_datablock<std::vector<double>>(
               sample, "cluster_abundance", "centered_cluster_biased_counts"),
             centered_count_errors = get_datablock<std::vector<double>>(
               sample, "cluster_abundance", "centered_cluster_count_errors"),
             miscentered_count_errors = get_datablock<std::vector<double>>(
               sample, "cluster_abundance", "centered_cluster_count_errors");

  // Start with Poisson errors and integration errors
  for (auto i = 0u; i < ncounts; i++) {
    centered_cluster_covariance[(ncounts * i) + i] = centered_counts[i];
    miscentered_cluster_covariance[(ncounts * i) + i] = miscentered_counts[i];
    centered_cluster_covariance[(ncounts * i) + i] +=
      centered_count_errors[i] * centered_count_errors[i];
    miscentered_cluster_covariance[(ncounts * i) + i] +=
      miscentered_count_errors[i] * miscentered_count_errors[i];
  }

  // Compute Sample Variance Covariance
  std::cout << "About to compute sigma^2...\n";
  const auto sigma_sq = sv.compute();
  std::cout << "Just computed sigma^2!\n";

  // Print results
  for (auto const& row : sigma_sq) {
    for (auto const& item : row) std::cerr << item << " ";
    std::cerr << '\n';
  }

  // Apply Sample Variance
  // C_{ij} = Nb_i * Nb_j * \sigma_{ij}^2
  for (auto i = 0u; i < ncounts; i++) {
    for (auto j = 0u; j < ncounts; j++) {
      const auto redshift_i = i % NREDSHIFT;
      const auto redshift_j = j % NREDSHIFT;
      const auto sample_variance_cen = centered_biased_counts[i] *
                                       centered_biased_counts[j] *
                                       sigma_sq[redshift_i][redshift_j];
      const auto sample_variance_miscen = miscentered_biased_counts[i] *
                                          miscentered_biased_counts[j] *
                                          sigma_sq[redshift_i][redshift_j];
      centered_cluster_covariance[(ncounts * i) + j] += sample_variance_cen;
      miscentered_cluster_covariance[(ncounts * i) + j] +=
        sample_variance_miscen;
    }
  }

  // Store Covariance
  std::vector<std::size_t> extents{{ncounts, ncounts}};
  sample.put_val<cosmosis::ndarray<double>>(
    "cluster_abundance",
    "centered_cluster_count_covariance",
    {centered_cluster_covariance, extents});
  sample.put_val<cosmosis::ndarray<double>>(
    "cluster_abundance",
    "miscentered_cluster_count_covariance",
    {miscentered_cluster_covariance, extents});

  return DBS_SUCCESS;
}

// Module interface implementation

extern "C" {

void*
setup(cosmosis::DataBlock* options)
{
  return new ClusterCovarianceModule(*options);
}

DATABLOCK_STATUS
execute(cosmosis::DataBlock* block, void* config)
{
  auto mod = static_cast<ClusterCovarianceModule*>(config);

  return mod->execute<DefaultModels::OMEGA_Z>(*block);
}

int
cleanup(void* config)
{
  delete static_cast<ClusterCovarianceModule*>(config);
  return 0;
}
}
