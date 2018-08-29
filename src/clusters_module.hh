#ifndef Y3CLUSTER_CLUSTERS_MODULE_HH
#define Y3CLUSTER_CLUSTERS_MODULE_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "/cosmosis/cubacpp/cubacpp/cubacpp.hh"
#include "gamma_t.hh"

namespace y3_cluster {
  template <class MODELS, std::size_t NRADII, std::size_t NRICHNESS = 1, std::size_t NREDSHIFT = 1>
  class ClustersModule {
    std::array<double, NRADII> radii_bins;
    std::array<y3_cluster::IntegrationRange, NRICHNESS> lo_bins;
    std::array<y3_cluster::IntegrationRange, NREDSHIFT> zo_bins;

  public:
    explicit ClustersModule(cosmosis::DataBlock& config);
    void execute(cosmosis::DataBlock& sample);
  };
}

template <std::size_t NRanges>
static inline std::array<y3_cluster::IntegrationRange, NRanges>
_get_ranges(cosmosis::DataBlock db, std::string name)
{
    std::array<std::size_t, NRanges> r;
    for (auto i = 0u; i < NRanges; i++)
        r[i] = i;

    auto edges = get_datablock<std::vector<double>>(db, OPTION_SECTION,
                                                    (name + "_bins").c_str());

    // TODO:
    //    This is not ideal. In the future we expect the number of bins
    //    to be configurable at runtime, not a template parameter.
    if (edges.size() != (NRanges + 1)) {
        throw std::runtime_error("Wrong number of edges for bins: " + name
                                 + " (expected " + std::to_string(NRanges + 1) + ")");
    }

    return y3_cluster::transform(r, [&](std::size_t i) {
              return y3_cluster::IntegrationRange{edges[i], edges[i + 1]};
            });
}

template <class MODELS, std::size_t NRADII, std::size_t NRICHNESS, std::size_t NREDSHIFT>
y3_cluster::ClustersModule<MODELS, NRADII, NRICHNESS, NREDSHIFT>::ClustersModule(cosmosis::DataBlock& config)
  // TODO: Possibly set up any optional parameters, like integration params?
  : lo_bins(_get_ranges<NRICHNESS>(config, "lo"))
  , zo_bins(_get_ranges<NREDSHIFT>(config, "zo"))
{
    auto radii = get_datablock<std::vector<double>>(config, OPTION_SECTION, "radii_bins");
    if (radii.size() != NRADII)
        throw std::runtime_error("Wrong number of radii! (expected " + std::to_string(NRADII) + ")");
    for (auto i = 0u; i < NRADII; i++)
        radii_bins[i] = radii[0];
}

// TODO:
// - Combine centered and miscentered terms
// - Compute covariance matrix
// - Compute likelihood - done in Michel's module
//
// Then you've got science!
template <class MODELS, std::size_t NRADII, std::size_t NRICHNESS, std::size_t NREDSHIFT>
void
y3_cluster::ClustersModule<MODELS, NRADII, NRICHNESS, NREDSHIFT>::execute(cosmosis::DataBlock& sample)
{
  // TODO should these be configurable from CosmoSIS?
  double const epsrel = 1.0e-3;
  double const epsabs = 1.0e-12;
  cubacpp::Cuhre c;
  c.maxeval = 100000000;

  // Initialize our big computation stuff up here, so any DataBlock access
  // errors happen up front
  std::vector<y3_cluster::IntegrationRange> zo_bins_vec;
  for (auto i = 0u; i < zo_bins.size(); i++)
    zo_bins_vec.push_back(zo_bins[i]);
  typename MODELS::OMEGA_Z omega_z(sample);
  typename MODELS::SAMPLE_VARIANCE sv(sample, omega_z, zo_bins_vec);
  Gamma_T_Integrand<MODELS, NRADII, NRICHNESS, NREDSHIFT> integrand(sample, radii_bins, lo_bins, zo_bins);

  // Create cluster count covariance - initialized to zeroes
  std::size_t ncounts = NRICHNESS * NREDSHIFT;
  std::vector<double> centered_cluster_count_covariance(ncounts * ncounts),
                      miscentered_cluster_count_covariance(ncounts * ncounts);
  for (auto i = 0u; i < ncounts*ncounts; i++) {
    centered_cluster_count_covariance[i] = 0.0;
    miscentered_cluster_count_covariance[i] = 0.0;
  }

  // Compute abundance counts and gamma_t's
  auto [centered_result, binned_centered_result] = integrand.integrate_centered(c, epsrel, epsabs);
  if (centered_result.status != 0)
    std::cerr << "WARNING: Centered result did not converge!\n";

  auto [miscentered_result, binned_miscentered_result] = integrand.integrate_miscentered(c, epsrel, epsabs);
  if (miscentered_result.status != 0)
    std::cerr << "WARNING: Miscentered result did not converge!\n";

  std::cout << "Centered:\n" << centered_result;
  std::cout << "Miscentered:\n" << miscentered_result;

  std::vector<double> centered_gamma_ts,
                      centered_cluster_counts,
                      miscentered_gamma_ts,
                      miscentered_cluster_counts;

  // Sort abundance counts and gamma_t's
  for (auto i = 0u; i < binned_centered_result.size(); i++) {
      const auto& bin = binned_centered_result[i];
      centered_cluster_counts.push_back(bin.N);
      // Store Poisson variance
      centered_cluster_count_covariance[(ncounts * i) + i] = bin.N;
      // Covariance due to integrator error
      centered_cluster_count_covariance[(ncounts * i) + i] += bin.N_error * bin.N_error;
      for (auto i = 0u; i < NRADII; i++)
          centered_gamma_ts.push_back(bin.gamma_ts[i]);
  }

  for (auto i = 0u; i < binned_miscentered_result.size(); i++) {
      const auto& bin = binned_miscentered_result[i];
      miscentered_cluster_counts.push_back(bin.N);
      // Store Poisson variance
      miscentered_cluster_count_covariance[(ncounts * i) + i] = bin.N;
      // Covariance due to integrator error
      miscentered_cluster_count_covariance[(ncounts * i) + i] += bin.N_error * bin.N_error;
      for (auto i = 0u; i < NRADII; i++)
          miscentered_gamma_ts.push_back(bin.gamma_ts[i]);
  }

  // Store abundance counts and gamma_t's
  sample.put_val<std::vector<double>>("cluster_abundance", "centered_gamma_ts",
                                      centered_gamma_ts);
  sample.put_val<std::vector<double>>("cluster_abundance", "centered_cluster_counts",
                                      centered_cluster_counts);
  sample.put_val<std::vector<double>>("cluster_abundance", "miscentered_gamma_ts",
                                      miscentered_gamma_ts);
  sample.put_val<std::vector<double>>("cluster_abundance", "miscentered_cluster_counts",
                                      miscentered_cluster_counts);

  // Compute Sample Variance Covariance
  std::cout << "About to compute sigma^2...\n";
  const auto sigma_sq = sv.compute();
  std::cout << "Just computed sigma^2!\n";

  for (const auto row : sigma_sq) {
    for (const auto item : row)
      std::cerr << item << " ";
    std::cerr << '\n';
  }

  // Apply Sample Variance
  // C_{ij} = Nb_i * Nb_j * \sigma_{ij}^2
  for (auto i = 0u; i < ncounts; i++) {
    for (auto j = 0u; j < ncounts; j++) {
      const auto redshift_i = i % NREDSHIFT;
      const auto redshift_j = j % NREDSHIFT;
      const auto sample_variance_cen = binned_centered_result[i].Nb * binned_centered_result[j].Nb
                                     * sigma_sq[redshift_i][redshift_j];
      const auto sample_variance_miscen = binned_miscentered_result[i].Nb * binned_miscentered_result[j].Nb
                                        * sigma_sq[redshift_i][redshift_j];
      centered_cluster_count_covariance[(ncounts * i) + j] += sample_variance_cen;
      miscentered_cluster_count_covariance[(ncounts * i) + j] += sample_variance_miscen;
    }
  }

  // Store Covariance
  std::vector<std::size_t> extents{{ncounts, ncounts}};
  sample.put_val<cosmosis::ndarray<double>>("cluster_abundance", "centered_cluster_count_covariance",
                                            {centered_cluster_count_covariance, extents});
  sample.put_val<cosmosis::ndarray<double>>("cluster_abundance", "miscentered_cluster_count_covariance",
                                            {miscentered_cluster_count_covariance, extents});

  // TODO gamma_t covariance
}

#endif
