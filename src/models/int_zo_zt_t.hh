#ifndef Y3_CLUSTER_INT_ZO_ZT_HH
#define Y3_CLUSTER_INT_ZO_ZT_HH

#include "cosmosis/datablock/datablock.hh"
#include <cmath>

namespace y3_cluster {

  class INT_ZO_ZT_t {
  public:
    explicit INT_ZO_ZT_t(double sigma) : _sigma(sigma) {}

    explicit INT_ZO_ZT_t(cosmosis::DataBlock& sample)
      : _sigma(sample.view<double>("cluster_abundance", "zo_zt_sigma"))
    {}

    double
    operator()(double zomin, double zomax, double zt) const
    {
      // P(zo | zt) := (1.0 / sqrt(2pi) / sigma) * exp(- (zo - zt) * (zo - zt) /
      // (2 * sigma * sigma))
      //    (i.e., a standard gaussian)
      // So,
      //    \int P(zo|zt) d(zo), zo in [zomin, zomax]
      //     == (erf((zomax - zt) / base) - erf((zomin - zt) / base)) / 2
      using std::erf;
      double const sigma2 = 0.02129638 - 0.25085154 * zt + 1.11756659 * zt * zt -
                1.22925508 * zt * zt * zt;
      double const base = std::sqrt(2.0) * sigma2;
      return (erf((zomax - zt) / base) - erf((zomin - zt) / base)) / 2.0;
    }

  private:
    double _sigma;
  };
}

#endif
