#ifndef Y3_CLUSTER_PZSOURCE_GAUSSIAN_HH
#define Y3_CLUSTER_PZSOURCE_GAUSSIAN_HH

#include "cosmosis/datablock/datablock.hh"
#include "utils/primitives.hh"

namespace y3_cluster {
  class PZSOURCE_GAUSSIAN_t {
  public:
    // TODO take in pzsource distributions
    explicit PZSOURCE_GAUSSIAN_t(double zbin, double sigma)
      : _zbin(zbin), _sigma(sigma)
    {}

    explicit PZSOURCE_GAUSSIAN_t(cosmosis::DataBlock& sample)
      : _zbin(sample.view<double>("cluster_abundance", "pzsource_zbin"))
      , _sigma(sample.view<double>("cluster_abundance", "pzsource_sigma"))
    {}

    double
    operator()(std::size_t /* bin number */, double zs) const
    {
      return y3_cluster::gaussian(zs, _zbin, _sigma);
    }

    constexpr std::size_t
    nsources() const
    {
      return 1;
    }

  private:
    double _zbin;
    double _sigma;
  };
}

#endif
