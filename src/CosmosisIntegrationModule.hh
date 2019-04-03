#ifndef Y3_CLUSTER_COSMOSIS_INTEGRATION_MODULE_HH
#define Y3_CLUSTER_COSMOSIS_INTEGRATION_MODULE_HH

namespace y3_cluster
{
  template <class I>
  class CosmoSISIntegrationModule
  {
  public:
    using IntegrandType = I;
    
    void execute(DataBlock& sample) {
      igrand_.set_sample(sample);
      // for each mass bin
      // for each lambda bin
      auto results = integrate(igrand_); // for vector valued integrand
      
      // shove results that into sample.
    };
    
  private:
    
    IntegrandType igrand_;
  };
}

#endif