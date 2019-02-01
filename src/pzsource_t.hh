#ifndef Y3_CLUSTER_PZSOURCE_HH
#define Y3_CLUSTER_PZSOURCE_HH

#include <datablock_reader.hh>
#include <cmath>
#include "primitives.hh"

namespace y3_cluster {
  class PZSOURCE_t {
  public:
    // TODO take in pzsource distributions
    explicit PZSOURCE_t(std::vector<std::shared_ptr<y3_cluster::Interp1D const>> pzsources)
        : pzsources(pzsources)
    {}

    explicit PZSOURCE_t(cosmosis::DataBlock&,
                        std::vector<std::shared_ptr<y3_cluster::Interp1D const>> pzsources)
      : pzsources(pzsources)
    {}

    double
    operator()(std::size_t bin_no, double zs) const
    {
      return pzsources[bin_no]->eval(zs);
    }

    std::size_t
    nsources() const
    {
      return pzsources.size();
    }

  private:
    // TODO move this to a different pzsource_t implementation?
    std::vector<std::shared_ptr<y3_cluster::Interp1D const>> pzsources;
  };
}

#endif
