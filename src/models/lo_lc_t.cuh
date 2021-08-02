#ifndef LO_LC_T_CUH
#define LO_LC_T_CUH

namespace quad {

  __device__ inline double
  gaussian(double x, double mu, double sigma)
  {
    double const z = (x - mu) / sigma;
    return exp(-z * z / 2.) * 0.3989422804014327 / sigma;
  }

}

class LO_LC_t {
public:
  LO_LC_t() = default;

  LO_LC_t(double alpha, double a, double b, double R_lambda)
    : _alpha(alpha), _a(a), _b(b), _R_lambda(R_lambda)
  {}

  __device__ double
  operator()(double lo, double lc, double R_mis) const
  {
    /* eq. (35) */
    double x = R_mis / _R_lambda;
    double y = lo / lc;
    double mu_y = exp(-x * x / _alpha / _alpha);
    double sigma_y = _a * atan(_b * x);
    // Need 1/lc scaling for total probability = 1
    return quad::gaussian(y, mu_y, sigma_y) / lc;
  }

private:
  double _alpha = 0.0;
  double _a = 0.0;
  double _b = 0.0;
  double _R_lambda = 0.0;
};

#endif
