#ifndef Y3_CLUSTER_CPP_MOR_T_HH
#define Y3_CLUSTER_CPP_MOR_T_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "test/mz_power_law.hh"

#include <cmath>

namespace y3_cluster {

  class MOR_t {
  public:
    MOR_t(mz_power_law lambda, double sigma, double alpha)
      : _lambda(lambda), _sigma(sigma), _alpha(alpha)
    {}

    explicit MOR_t(cosmosis::DataBlock& sample)
      : _lambda([](cosmosis::DataBlock& x) {
        double A, B, C;
        x.get_val<double>("MOR_params", "A", A);
        x.get_val<double>("MOR_params", "B", B);
        x.get_val<double>("MOR_params", "C", C);
        return mz_power_law{A, B, C};
      }(sample))
    {
      sample.get_val<double>("MOR_params", "sigma", _sigma);
      sample.get_val<double>("MOR_params", "alpha", _alpha);
    }

    double
    operator()(double lt, double lnM, double zt) const
    {
      /* eq. (34) */
      double const ltm = _lambda(lnM, zt);
      double const x = lt - ltm;
      double const erfarg = -1.0 * _alpha * (x) / (std::sqrt(2.) * _sigma);
      double const erfterm = std::erfc(erfarg);
      return y3_cluster::gaussian(x, 0.0, _sigma) * erfterm;
    }

  private:
    mz_power_law _lambda;
    double _sigma;
    double _alpha;
  };
}

#endif
