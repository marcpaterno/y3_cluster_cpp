#ifndef Y3CLUSTER_CLUSTERS_MODULE_HH
#define Y3CLUSTER_CLUSTERS_MODULE_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "/cosmosis/cubacpp/cubacpp/cubacpp.hh"
#include "gamma_t.hh"
#include <datablock_reader.hh>
#include "gperftools/profiler.h"

namespace y3_cluster {
  template <class MODELS>
  class ClustersModule {
    bool clusters_module_on;
    ProfilerOptions options;
    int verbosity = 1;
    bool profile = false;
    std::vector<double> radii_bins;
    std::vector<y3_cluster::IntegrationRange> lo_bins;
    std::vector<y3_cluster::IntegrationRange> zo_bins;

  public:
    explicit ClustersModule(cosmosis::DataBlock& config);
    ~ClustersModule() {
      ProfilerFlush();
      ProfilerStop();
    };
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

extern "C" {
  int filter_in_thread(void *arg) {
    volatile bool *on = reinterpret_cast<bool *>(arg);
    return (*on) ? 1 : 0;
  }
}

template <class MODELS>
y3_cluster::ClustersModule<MODELS>::ClustersModule(cosmosis::DataBlock& config)
  // TODO: Possibly set up any optional parameters, like integration params?
  : clusters_module_on(false)
  , options{filter_in_thread, &clusters_module_on}
  , verbosity(get_datablock<int>(config, "cluster_abundance", "verbosity", 1))
  , profile(get_datablock<bool>(config, "cluster_abundance", "profile", false))
  , radii_bins(get_datablock<std::vector<double>>(config, OPTION_SECTION, "radii_bins"))
  , lo_bins(_get_ranges(config, "lo"))
  , zo_bins(_get_ranges(config, "zo"))
{
  if (profile)
    ProfilerStartWithOptions("/cosmosis/cosmosis-standard-library/y3_cluster_cpp/cosmosis_run_dump.txt",
                             &options);
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
  // FIXME: Just a test placeholder! Should these come from CosmoSIS?
  double const epsrel = 1.0e-3;
  double const epsabs = 1.0e-12;
  cubacores(0,0);
  cubacpp::Cuhre c;
  c.maxeval = 100000000;

  // Initialize our big computation stuff up here, so any DataBlock access
  // errors happen up front
  auto integrand = Gamma_T_Integrand<MODELS>
                     {sample, radii_bins, lo_bins, zo_bins};

  // Compute abundance counts and gamma_t's
  // Profile if necessary
  clusters_module_on = true;
  auto centered_result = integrand.integrate_centered(c, epsrel, epsabs);
  auto miscentered_result = integrand.integrate_miscentered(c, epsrel, epsabs);
  clusters_module_on = false;

  if (centered_result.status != 0) {
    std::cout << "h0 " << get_datablock<double>(sample, "cosmological_parameters", "h0") << "\n";
    std::cout << "Omega_m " << get_datablock<double>(sample, "cosmological_parameters", "omega_m") << "\n";
    std::cout << "Omega_b " << get_datablock<double>(sample, "cosmological_parameters", "omega_b") << "\n";
    std::cout << "Omega_nu " << get_datablock<double>(sample, "cosmological_parameters", "omega_nu") << "\n";
    std::cout << "log1e10As " << get_datablock<double>(sample, "cosmological_parameters", "log1e10As") << "\n";
    std::cout << "n_s " << get_datablock<double>(sample, "cosmological_parameters", "n_s") << "\n";
    std::cout << "mor_A " << get_datablock<double>(sample, "cluster_abundance", "mor_A") << "\n";
    std::cout << "mor_B " << get_datablock<double>(sample, "cluster_abundance", "mor_B") << "\n";
    std::cout << "mor_alpha " << get_datablock<double>(sample, "cluster_abundance", "mor_alpha") << "\n";
    std::cout << "mor_sigma " << get_datablock<double>(sample, "cluster_abundance", "mor_sigma") << "\n";
    std::cout << "hmf_s " << get_datablock<double>(sample, "cluster_abundance", "hmf_s") << "\n";
    std::cout << "hmf_q " << get_datablock<double>(sample, "cluster_abundance", "hmf_q") << "\n";
    throw std::runtime_error("Centered integration did not converge!");
  }
  if (miscentered_result.status != 0) {
    std::cout << "h0 " << get_datablock<double>(sample, "cosmological_parameters", "h0") << "\n";
    std::cout << "Omega_m " << get_datablock<double>(sample, "cosmological_parameters", "omega_m") << "\n";
    std::cout << "Omega_b " << get_datablock<double>(sample, "cosmological_parameters", "omega_b") << "\n";
    std::cout << "Omega_nu " << get_datablock<double>(sample, "cosmological_parameters", "omega_nu") << "\n";
    std::cout << "log1e10As " << get_datablock<double>(sample, "cosmological_parameters", "log1e10As") << "\n";
    std::cout << "n_s " << get_datablock<double>(sample, "cosmological_parameters", "n_s") << "\n";
    std::cout << "mor_A " << get_datablock<double>(sample, "cluster_abundance", "mor_A") << "\n";
    std::cout << "mor_B " << get_datablock<double>(sample, "cluster_abundance", "mor_B") << "\n";
    std::cout << "mor_alpha " << get_datablock<double>(sample, "cluster_abundance", "mor_alpha") << "\n";
    std::cout << "mor_sigma " << get_datablock<double>(sample, "cluster_abundance", "mor_sigma") << "\n";
    std::cout << "hmf_s " << get_datablock<double>(sample, "cluster_abundance", "hmf_s") << "\n";
    std::cout << "hmf_q " << get_datablock<double>(sample, "cluster_abundance", "hmf_q") << "\n";
    throw std::runtime_error("Miscentered integration did not converge!");
  }
  if (verbosity > 0 ) {
    std::cout << "Centered:\n" << centered_result;
    std::cout << "Miscentered:\n" << miscentered_result;
  } else {
    std::cout << "Centered and miscentered integration\n";
  }

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
    for (const auto gt : bin.gamma_ts)
      centered_gamma_ts.push_back(gt);
  }

  // Sort miscentered abundance counts/biased counts/error bars and gamma_t's
  for (const auto& bin : miscentered_result) {
    miscentered_cluster_counts.push_back(bin.N);
    miscentered_count_errors.push_back(bin.N_error);
    miscentered_biased_counts.push_back(bin.Nb);
    for (const auto gt : bin.gamma_ts)
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
