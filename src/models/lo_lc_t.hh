#ifndef Y3_CLUSTER_CPP_LO_LC_T_HH
#define Y3_CLUSTER_CPP_LO_LC_T_HH

#include "utils/datablock_reader.hh"
#include "utils/primitives.hh"

#include <cmath>

// In Yuanyaun's section 9 summary of the Lighthouse document
// https://www.overleaf.com/project/5c378b07f882d02f5b8c90e2
// this is equation 53, labeled P(\lambda_o|\lambda_c,R_cent)

namespace y3_cluster {

  class LO_LC_t {
  public:
    LO_LC_t(double alpha, double a, double b, double R_lambda)
      : _alpha(alpha), _a(a), _b(b), _R_lambda(R_lambda)
    {}

    explicit LO_LC_t(cosmosis::DataBlock& sample)
      : _alpha(sample.view<double>("cluster_abundance", "LO_LC_alpha"))
      , _a(sample.view<double>("cluster_abundance", "LO_LC_a"))
      , _b(sample.view<double>("cluster_abundance", "LO_LC_b"))
      , _R_lambda(sample.view<double>("cluster_abundance", "LO_LC_R_lambda"))
    {}

    double
    operator()(double lo, double lc, double R_mis) const
    {
      //  https://www.overleaf.com/project/5c378b07f882d02f5b8c90e2 
      // equation 35, and  equation 53 in Yuanyuan's clean up section 
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
