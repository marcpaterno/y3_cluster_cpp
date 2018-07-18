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
  : lo_bins(_get_ranges<NREDSHIFT>(config, "lo"))
  , zo_bins(_get_ranges<NREDSHIFT>(config, "zo"))
{
    auto radii = get_datablock<std::vector<double>>(config, OPTION_SECTION, "radii_bins");
    if (radii.size() != NRADII)
        throw std::runtime_error("Wrong number of radii! (expected " + std::to_string(NRADII) + ")");
    for (auto i = 0u; i < NRADII; i++)
        radii_bins[i] = radii[0];
}

template <class MODELS, std::size_t NRADII, std::size_t NRICHNESS, std::size_t NREDSHIFT>
void
y3_cluster::ClustersModule<MODELS, NRADII, NRICHNESS, NREDSHIFT>::execute(cosmosis::DataBlock& sample)
{
  // FIXME: Just a test placeholder!
  // sample.put_val("cluster_abundance", "likelihood", val);
  double const epsrel = 1.0e-3;
  double const epsabs = 1.0e-12;
  cubacpp::Cuhre c;
  c.maxeval = 100000000;

  Gamma_T_Integrand<MODELS, NRADII, NRICHNESS, NREDSHIFT> integrand(sample, radii_bins, lo_bins, zo_bins);
  auto centered_result = integrand.integrate_centered(c, epsrel, epsabs);
  auto miscentered_result = integrand.integrate_miscentered(c, epsrel, epsabs);

  std::cout << "Centered:\n" << centered_result;
  std::cout << "Miscentered:\n" << miscentered_result;

  // TODO:
  // - Combine centered and miscentered terms
  // - Compute covariance matrix
  // - Compute likelihood
  //
  // Then you've got science!

  //sample.put_val("cluster_abundance", "likelihood", result);

//  time_integration(c, gti, epsrel, epsabs, "cuhre");
//  time_integration(ALG alg,
//  		 F f,
//  		 double epsrel,
//  		 double epsabs,
//  		 char const* algname)
//  {
//    auto start = std::chrono::high_resolution_clock::now();
//    auto res = alg.integrate(f, epsrel, epsabs);
//    auto stop = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> diff = stop - start;
//    std::cout << algname << ": " << res << " (" << diff.count() << "s)\n";
//  }
}

#endif
