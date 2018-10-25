#ifndef Y3_CLUSTER_ZO_ZT_HH
#define Y3_CLUSTER_ZO_ZT_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"

#include <datablock_reader.hh>
#include <cmath>

namespace y3_cluster {

  class ZO_ZT_t {
  public:
    explicit ZO_ZT_t(double sigma) : _sigma(sigma) {}

    explicit ZO_ZT_t(cosmosis::DataBlock& sample)
      : _sigma(get_datablock<double>(sample, "cluster_abundance", "zo_zt_sigma"))
    {}

    double
    operator()(double zomin, double zomax, double zt) const
    {
      // P(zo | zt) := (1.0 / sqrt(2pi) / sigma) * exp(- (zo - zt) * (zo - zt) / (2 * sigma * sigma))
      //    (i.e., a standard gaussian)
      // So,
      //    \int P(zo|zt) d(zo), zo in [zomin, zomax]
      //     == (erf((zomax - zt) / base) - erf((zomin - zt) / base)) / 2
      using std::erf;
      double base = std::sqrt(2) * _sigma;
      return (erf((zomax - zt) / base) -
              erf((zomin - zt) / base)) / 2.0;
    }

  private:
    double _sigma;
  };
}

#endif
