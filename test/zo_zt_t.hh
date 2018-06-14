#ifndef Y3_CLUSTER_ZO_ZT_HH
#define Y3_CLUSTER_ZO_ZT_HH

namespace y3_cluster {

  class ZO_ZT_t {
  public:
    explicit ZO_ZT_t(double sigma) : _sigma(sigma) {}

    explicit ZO_ZT_t(cosmosis::DataBlock& sample)
    {
      sample.get_val<double>("ZO_ZT_params", "sigma", _sigma);
    }

    double
    operator()(double zo, double zt) const
    {
      /* eq. (32) */
      double const x = zo - zt;
      return y3_cluster::gaussian(x, 0., _sigma);
    }

  private:
    double _sigma;
  };
}

#endif
