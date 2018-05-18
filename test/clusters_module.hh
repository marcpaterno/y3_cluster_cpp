#ifndef Y3CLUSTER_CLUSTERS_MODULE_HH
#define Y3CLUSTER_CLUSTERS_MODULE_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "gamma_t.hh"

namespace y3_cluster {
  template <class MODELS>
  class ClustersModule {
  public:
    explicit ClustersModule(cosmosis::DataBlock& config);
    void execute(cosmosis::DataBlock& sample) const;

  private:
    gamma_t<MODELS> integrand_;
  };
}

template <class MODELS>
y3_cluster::ClustersModule<MODELS>::ClustersModule(cosmosis::DataBlock& config)
  : integrand_(config)
{}

template <class MODELS>
void
y3_cluster::ClustersModule<MODELS>::execute(cosmosis::DataBlock& sample)
{
  auto val = integrand_.evaluate(sample);
  sample.put_val("x", "y", val);
}

#endif
