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

  private:
    Gamma_T_Integrand<MODELS, NRADII> _integrand;
  };
}

template <class MODELS, std::size_t NRADII>
y3_cluster::ClustersModule<MODELS, NRADII>::ClustersModule(cosmosis::DataBlock& config)
  : _integrand(config)
{}

template <class MODELS, std::size_t NRADII>
void
y3_cluster::ClustersModule<MODELS, NRADII>::execute(cosmosis::DataBlock& sample)
{
  // FIXME: evaluate isn't a thing, it should look like this from trivial_gamma_t.cc:
  //    auto res = alg.integrate(f, epsrel, epsabs);
  // auto val = integrand_.evaluate(sample);
  // Placeholder...
  // sample.put_val("gamma_t", "likelihood", val);
  // FIXME: Just a test placeholder!
  double const epsrel = 1.0e-3;
  double const epsabs = 1.0e-12;
  cubacpp::Cuhre c;
  c.maxeval = 100000000;

  auto centered_result = _integrand.integrate_centered(c, epsrel, epsabs);
  auto miscentered_result = _integrand.integrate_miscentered(c, epsrel, epsabs);

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
