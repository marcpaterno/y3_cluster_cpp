#ifndef Y3_CLUSTER_INTERP_1D_HH
#define Y3_CLUSTER_INTERP_1D_HH

#include "gsl/gsl_interp.h"
#include <array>
#include <cstddef>
#include <vector>

// Interp1D is used for linear interpolation in 1 dimension.
// It uses the GSL library to do the actual interpolation.
// Interp1D object allow extrapolation as well as supporting
// interpolation; no warnings or errors are given when
// extrapolating.
//
namespace y3_cluster {
  class Interp1D {
  public:
    template <std::size_t N>
    Interp1D(std::array<double, N> const& xs, std::array<double, N> const& ys);

  private:
    std::vector<double> xs_;
    std::vector<double> ys_;
  };
}

template <std::size_t N>
y3_cluster::Interp1D::Interp1D(std::array<double, N> const& xs,
                               std::array<double, N> const& ys)
  : xs_(begin(xs), end(xs)), ys_(begin(ys), end(ys))
{}
#endif
