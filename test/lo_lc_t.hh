#ifndef Y3_CLUSTER_CPP_LO_LC_T_HH
#define Y3_CLUSTER_CPP_LO_LC_T_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "test/primitives.hh"
#include <cmath>

namespace y3_cluster {

  class LO_LC_t {
  public:
    LO_LC_t(double alpha, double a, double b, double R_lambda)
      : _alpha(alpha), _a(a), _b(b), _R_lambda(R_lambda)
    {}

    explicit LO_LC_t(cosmosis::DataBlock& sample)
    {
      sample.get_val<double>("LO_LC_params", "alpha", _alpha);
      sample.get_val<double>("LO_LC_params", "a", _a);
      sample.get_val<double>("LO_LC_params", "b", _b);
      sample.get_val<double>("LO_LC_params", "R_lambda", _R_lambda);
    }

    double
    operator()(double lo, double lc, double R_mis) const
    {
      /* eq. (35) */
      double x = R_mis / _R_lambda;
      double y = lo / lc;
      double mu_y = std::exp(-x * x / _alpha / _alpha);
      double sigma_y = _a * std::atan(_b * x);
      // Need 1/lc scaling for total probability = 1
      return y3_cluster::gaussian(y, mu_y, sigma_y) / lc;
    }

  private:
    double _alpha;
    double _a;
    double _b;
    double _R_lambda;
  };
}

#endif
