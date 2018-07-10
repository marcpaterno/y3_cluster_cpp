#ifndef Y3CLUSTER_CLUSTERS_MODULE_HH
#define Y3CLUSTER_CLUSTERS_MODULE_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "/cosmosis/cubacpp/cubacpp/cubacpp.hh"
#include "gamma_t.hh"

namespace y3_cluster {
  template <class MODELS>
  class ClustersModule {
  public:
    explicit ClustersModule(cosmosis::DataBlock& config);
    void execute(cosmosis::DataBlock& sample) const;

  private:
    Gamma_T_Integrand<MODELS> _integrand;
  };
}

template <class MODELS>
y3_cluster::ClustersModule<MODELS>::ClustersModule(cosmosis::DataBlock& config)
  : _integrand(config)
{}

template <class MODELS>
void
y3_cluster::ClustersModule<MODELS>::execute(cosmosis::DataBlock& sample) const
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

  auto result = c.integrate(_integrand, epsrel, epsabs);

  std::cout << result;

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
