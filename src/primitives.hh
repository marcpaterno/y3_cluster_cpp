#ifndef Y3_CLUSTER_CPP_PRIMITIVES_HH
#define Y3_CLUSTER_CPP_PRIMITIVES_HH

#include <polynomial.hh>

#include <cmath>
#include <fstream>
#include <iostream>
// primitives.hh contains a few commonly-used mathematical primitives.

namespace y3_cluster {

  double constexpr pi() { return 4. * std::atan(1.0); };

  // TODO: Get higher precision!
  double constexpr c() { return 9.71561e-15;}; // Mpc/s

  // TODO: Get higher precision!
  double constexpr G() { return 4.517e-48;}; // Mpc^3 / M_sol / s^2

  double constexpr invsqrt2pi() { return 1. / std::sqrt(2. * pi()); };

  inline double
  gaussian(double x, double mu, double sigma)
  {
    double const z = (x - mu) / sigma;
    return std::exp(-z * z / 2.) * 0.3989422804014327 / sigma;
  }

  namespace {
    // Tail recursive helper for `integer_pow`
    constexpr double
    do_integer_pow(const double accumulator, const double n, const unsigned pow)
    {
      if (pow == 0)
        return accumulator;
      if ((pow % 2) == 0)
        return do_integer_pow(accumulator, n * n, pow / 2);
      return do_integer_pow(accumulator * n, n, pow - 1);
    }
  }

  // In C++ >= 11, the std::pow does not optimize for integers :/
  constexpr double
  integer_pow(double n, int pow)
  {
    if (pow == 0)
      return 1;
    if (pow < 0)
      return do_integer_pow(1, 1.0 / n, 0 - pow);
    return do_integer_pow(1, n, pow);
  }

  /*
  double constexpr
  erf_approx(const double x)
  {
	  // Error function approximation obtained here:
	  // https://en.wikipedia.org/wiki/Error_function#Approximation_with_elementary_functions
	  // Valid to within 3e-7
      constexpr polynomial<7> poly{{0.0000430638,
                                    0.0002765672,
                                    0.0001520143,
                                    0.0092705272,
                                    0.0422820123,
                                    0.0705230784,
                                    1}};

      double sign = x < 0 ? -1 : 1,
             res = poly(sign * x);

      // Sixteenth power => (2^4) power
      res *= res;
      res *= res;
      res *= res;
      res *= res;

      return sign * (1 - 1 / res);
  }
  */
}

#endif
