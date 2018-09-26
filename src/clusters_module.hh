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
    auto const  edges  =  get_datablock<std::vector<double>>(db, OPTION_SECTION,
                                                             (name + "_bins").c_str());

    auto  ret  =  std::vector<y3_cluster::IntegrationRange> (edges.size () - 1);

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
    auto  const  radii  =  get_datablock<std::vector<double>>(config, OPTION_SECTION, "radii_bins");

    radii_bins  =  std::vector<double>  (radii.size (),  radii [0]);
}

template <class MODELS>
void
y3_cluster::ClustersModule<MODELS>::execute(cosmosis::DataBlock& sample)
{
  // FIXME: Just a test placeholder!
  // sample.put_val("cluster_abundance", "likelihood", val);
  double const epsrel = 1.0e-3;
  double const epsabs = 1.0e-12;
  cubacpp::Cuhre c;
  c.maxeval = 100000000;

  auto integrand = Gamma_T_Integrand<MODELS>
                     {sample, radii_bins, lo_bins, zo_bins};

  auto centered_result = integrand.integrate_centered(c, epsrel, epsabs);
  auto miscentered_result = integrand.integrate_miscentered(c, epsrel, epsabs);

  if (centered_result.status != 0)
    throw std::runtime_error("Centered integration did not converge!");
  if (miscentered_result.status != 0)
    throw std::runtime_error("Centered integration did not converge!");

  std::cout << "Centered:\n" << centered_result;
  std::cout << "Miscentered:\n" << miscentered_result;

  std::vector<double> centered_gamma_ts,
                      centered_cluster_counts,
                      miscentered_gamma_ts,
                      miscentered_cluster_counts;

  for (const auto& bin : centered_result) {
      centered_cluster_counts.push_back(bin.N);
      //centered_cluster_counts.push_back(bin.Nw);
      for (auto i = 0u; i < radii_bins.size (); i++)
          centered_gamma_ts.push_back(bin.gamma_ts[i]);
  }

  for (const auto& bin : miscentered_result) {
      miscentered_cluster_counts.push_back(bin.N);
      //miscentered_cluster_counts.push_back(bin.Nw);
      for (auto i = 0u; i < radii_bins.size (); i++)
          miscentered_gamma_ts.push_back(bin.gamma_ts[i]);
  }

  sample.put_val<std::vector<double>>("cluster_abundance", "centered_gamma_ts",
                                      centered_gamma_ts);
  sample.put_val<std::vector<double>>("cluster_abundance", "centered_cluster_counts",
                                      centered_cluster_counts);
  sample.put_val<std::vector<double>>("cluster_abundance", "miscentered_gamma_ts",
                                      miscentered_gamma_ts);
  sample.put_val<std::vector<double>>("cluster_abundance", "miscentered_cluster_counts",
                                      miscentered_cluster_counts);

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
