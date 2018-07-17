#ifndef Y3CLUSTER_CLUSTERS_MODULE_HH
#define Y3CLUSTER_CLUSTERS_MODULE_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "/cosmosis/cubacpp/cubacpp/cubacpp.hh"
#include "gamma_t.hh"

namespace y3_cluster {
  template <class MODELS, std::size_t NRADII>
  class ClustersModule {
  public:
    explicit ClustersModule(cosmosis::DataBlock& config);
    void execute(cosmosis::DataBlock& sample);
  };
}

template <class MODELS, std::size_t NRADII>
y3_cluster::ClustersModule<MODELS, NRADII>::ClustersModule(cosmosis::DataBlock& config)
  // TODO: Possibly set up any optional parameters, like integration params?
{}

template <class MODELS, std::size_t NRADII>
void
y3_cluster::ClustersModule<MODELS, NRADII>::execute(cosmosis::DataBlock& sample)
{
  // FIXME: Just a test placeholder!
  // sample.put_val("gamma_t", "likelihood", val);
  double const epsrel = 1.0e-3;
  double const epsabs = 1.0e-12;
  cubacpp::Cuhre c;
  c.maxeval = 100000000;

  Gamma_T_Integrand<MODELS, NRADII> integrand(sample);
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

  //sample.put_val("gamma_t", "likelihood", result);

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
