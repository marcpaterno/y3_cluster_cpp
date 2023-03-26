#ifndef Y3_CLUSTER_CPP_PRIMITIVES_CUH
#define Y3_CLUSTER_CPP_PRIMITIVES_CUH

// primitives.cuh contains a few commonly-used mathematical primitives.

#include <limits>

namespace y3_cuda {

  namespace detail {
    double constexpr sqrtNewtonRaphson(double x, double curr, double prev)
    {
      return curr == prev ? curr :
                            sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
    }
  } // namespace detail

  // Constexpr version of the square root.
  // For a finite and non-negative value of "x", returns an approximation for
  // the square root of "x" Otherwise, returns NaN. Adapted from
  // https://stackoverflow.com/questions/8622256/in-c11-is-sqrt-defined-as-constexpr
  double constexpr sqrt(double x)
  {
    return (x >= 0.0 && isfinite(x)) ? detail::sqrtNewtonRaphson(x, x, 0) :
                                       std::numeric_limits<double>::quiet_NaN();
  }
  // This is the double-precision hexidecimal floating point value
  // closes to pi.
  double constexpr pi()
  {
    return 0x1.5bf0a8b145769p+1;
  }

  // TODO: Get higher precision!
  // Source: astropy's constants and unit conversion
  double constexpr c()
  {
    return 9.71561189e-15;
  } // Mpc/s

  // TODO: Get higher precision!
  // Source: astropy's constants and unit conversion
  double constexpr G()
  {
    return 4.51710305e-48;
  } // Mpc^3 / M_sol / s^2

  double constexpr invsqrt2pi()
  {
    return 1. / sqrt(2. * pi());
  }

  inline __device__ __host__ double
  gaussian(double x, double mu, double sigma)
  {
    double const z = (x - mu) / sigma;
    return exp(-z * z / 2.) * 0.3989422804014327 / sigma;
  }
  
  
  // Bivariate Gaussian Implementation
  // https://online.stat.psu.edu/stat505/lesson/4/4.2
  inline __device__ __host__ double
  gaussian2d(double x1, double x2, double mu1, double mu2, double sigma1, double sigma2, double corr12)
  {
    double const z1 = (x1 - mu1) / sigma1;
    double const z2 = (x2 - mu2) / sigma2;
    double const z12 = z1*z2

    double const rhosqr = max(min(1. - corr12 * corr12, 0.), 1.)
    double const arg = -(z1 * z1 + z2 * z2 + 2 * corr12 * z12)/2./rhosqr
    double const sigma = sigma1 * sigma2 * sqrt(rhosqr)

    return exp(arg) * 0.3989422804014327 / sigma;
  }
  namespace {
    // Tail recursive helper for `integer_pow`
    constexpr double
    do_integer_pow(const double accumulator, const double n, const unsigned pow)
    {
      if (pow == 0) return accumulator;
      if ((pow % 2) == 0) return do_integer_pow(accumulator, n * n, pow / 2);
      return do_integer_pow(accumulator * n, n, pow - 1);
    }
  } // namespace

  // In C++ >= 11, the std::pow does not optimize for integers :/
  constexpr double
  integer_pow(double n, int pow)
  {
    if (pow == 0) return 1;
    if (pow < 0) return do_integer_pow(1, 1.0 / n, 0 - pow);
    return do_integer_pow(1, n, pow);
  }
} // namespace y3_cuda

#endif
