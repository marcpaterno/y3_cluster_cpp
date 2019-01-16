#ifndef Y3_CLUSTER_PZSOURCE_HH
#define Y3_CLUSTER_PZSOURCE_HH

#include <datablock_reader.hh>
#include <cmath>
#include "primitives.hh"

namespace y3_cluster {

  class PZSOURCE_t {
  public:
    explicit PZSOURCE_t(double zbin, double sigma) : _zbin(zbin), _sigma(sigma) {}

    explicit PZSOURCE_t(cosmosis::DataBlock& sample)
      : _zbin(get_datablock<double>(sample, "cluster_abundance", "pzsource_zbin")),
        _sigma(get_datablock<double>(sample, "cluster_abundance", "pzsource_sigma"))
    {}

    double
    operator()(double zs) const
    {
      return y3_cluster::gaussian(zs, _zbin, _sigma);
    }

  private:
    double _zbin;
    double _sigma;
  };
}

#endif
