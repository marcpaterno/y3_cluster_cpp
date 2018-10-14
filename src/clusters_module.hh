#ifndef Y3CLUSTER_CLUSTERS_MODULE_HH
#define Y3CLUSTER_CLUSTERS_MODULE_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "/cosmosis/cubacpp/cubacpp/cubacpp.hh"
#include "gamma_t.hh"

namespace y3_cluster {
  template <class MODELS>
  class ClustersModule {
    std::vector<double> radii_bins;
    std::vector<y3_cluster::IntegrationRange> lo_bins;
    std::vector<y3_cluster::IntegrationRange> zo_bins;

  public:
    explicit ClustersModule(cosmosis::DataBlock& config);
    void execute(cosmosis::DataBlock& sample);
  };
}


static inline std::vector<y3_cluster::IntegrationRange>
_get_ranges(cosmosis::DataBlock db, std::string name)
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

template <class MODELS>
y3_cluster::ClustersModule<MODELS>::ClustersModule(cosmosis::DataBlock& config)
  // TODO: Possibly set up any optional parameters, like integration params?
  : lo_bins(_get_ranges(config, "lo"))
  , zo_bins(_get_ranges(config, "zo"))
{
  auto const radii = get_datablock<std::vector<double>>(config, OPTION_SECTION, "radii_bins");

  radii_bins = std::vector<double>  (radii.size (),  radii [0]);
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
  std::size_t const NRICHNESS = lo_bins.size(),
                    NREDSHIFT = zo_bins.size(),
                    NRADII = radii_bins.size();

  // FIXME: Just a test placeholder! Should these come from CosmoSIS?
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
  auto integrand = Gamma_T_Integrand<MODELS>
                     {sample, radii_bins, lo_bins, zo_bins};

  // Create cluster count covariance - initialized to zeroes
  std::size_t ncounts = NRICHNESS * NREDSHIFT;
  std::vector<double> centered_cluster_count_covariance(ncounts * ncounts),
                      miscentered_cluster_count_covariance(ncounts * ncounts);
  for (auto i = 0u; i < ncounts*ncounts; i++) {
    centered_cluster_count_covariance[i] = 0.0;
    miscentered_cluster_count_covariance[i] = 0.0;
  }

  // Compute abundance counts and gamma_t's
  auto centered_result = integrand.integrate_centered(c, epsrel, epsabs);
  auto miscentered_result = integrand.integrate_miscentered(c, epsrel, epsabs);

  if (centered_result.status != 0)
    throw std::runtime_error("Centered integration did not converge!");
  if (miscentered_result.status != 0)
    throw std::runtime_error("Miscentered integration did not converge!");

  std::cout << "Centered:\n" << centered_result;
  std::cout << "Miscentered:\n" << miscentered_result;

  std::vector<double> centered_gamma_ts,
                      centered_cluster_counts,
                      miscentered_gamma_ts,
                      miscentered_cluster_counts;

  // Sort abundance counts and gamma_t's
  for (auto i = 0u; i < centered_result.size(); i++) {
    const auto& bin = centered_result[i];
    centered_cluster_counts.push_back(bin.N);
    // Store Poisson variance
    centered_cluster_count_covariance[(ncounts * i) + i] = bin.N;
    // Covariance due to integrator error
    centered_cluster_count_covariance[(ncounts * i) + i] += bin.N_error * bin.N_error;
    for (auto i = 0u; i < NRADII; i++)
      centered_gamma_ts.push_back(bin.gamma_ts[i]);
  }

  for (auto i = 0u; i < miscentered_result.size(); i++) {
    const auto& bin = miscentered_result[i];
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
      const auto sample_variance_cen = centered_result[i].Nb * centered_result[j].Nb
                                     * sigma_sq[redshift_i][redshift_j];
      const auto sample_variance_miscen = miscentered_result[i].Nb * miscentered_result[j].Nb
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
