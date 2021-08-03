#ifndef Y3_CLUSTER_CPP_PRIMITIVES_CUH
#define Y3_CLUSTER_CPP_PRIMITIVES_CUH

// primitives.cuh contains a few commonly-used mathematical primitives.

#include <limits>

namespace y3_cuda{

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
    return (x >= 0.0 && isfinite(x)) ?
             detail::sqrtNewtonRaphson(x, x, 0) :
             std::numeric_limits<double>::quiet_NaN();
  }
  // This is the double-precision hexidecimal floating point value
  // closes to pi.
  double constexpr pi() { return 0x1.5bf0a8b145769p+1; }

  // TODO: Get higher precision!
  // Source: astropy's constants and unit conversion
  double constexpr c() { return 9.71561189e-15; } // Mpc/s

  // TODO: Get higher precision!
  // Source: astropy's constants and unit conversion
  double constexpr G() { return 4.51710305e-48; } // Mpc^3 / M_sol / s^2

  double constexpr invsqrt2pi() { return 1. / sqrt(2. * pi()); }

  inline 
  __device__ double
  gaussian(double x, double mu, double sigma)
  {
    double const z = (x - mu) / sigma;
    return exp(-z * z / 2.) * 0.3989422804014327 / sigma;
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
