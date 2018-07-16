#ifndef Y3_CLUSTER_CPP_PRIMITIVES_HH
#define Y3_CLUSTER_CPP_PRIMITIVES_HH

#include <cmath>
#include <fstream>
#include <iostream>
// primitives.hh contains a few commonly-used mathematical primitives.

namespace y3_cluster {

  double constexpr pi() { return 4. * std::atan(1.0); };

  double constexpr invsqrt2pi() { return 1. / std::sqrt(2. * pi()); };

  inline double
  gaussian(double x, double mu, double sigma)
  {
    double const z = (x - mu) / sigma;
    return std::exp(-z * z / 2.) * 0.3989422804014327 / sigma;
  }
}

#endif
