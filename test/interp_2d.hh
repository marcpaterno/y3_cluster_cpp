#ifndef Y3_CLUSTER_INTERP_2D_HH
#define Y3_CLUSTER_INTERP_2D_HH

#include "test/point_3d.hh"

#include "gsl/gsl_interp2d.h"
#include <array>
#include <cstddef>
#include <vector>

// Interp2D is used for linear interpolation in 1 dimension.
// It uses the GSL library to do the actual interpolation.
// Interp2D object allow extrapolation as well as supporting
// interpolation; no warnings or errors are given when
// extrapolating.
//
namespace y3_cluster {

  class Interp2D {
  public:
    // Interpolator created from two arrays; compiler assures they are of the
    // same length.
    template <std::size_t M, std::size_t N>
    Interp2D(std::array<double, M> const& xs,
             std::array<double, N> const& ys,
             std::array<std::array<double, M>, N> const& zs);

    // Interpolator created from vector of triplets: throws std::logic_error if
    // the points do not lie on a grid in (x,y) space, or if any values are NaN
    // or infinities. Any denormalized x- or y-values are flushed to zero. We
    // take the vector by value because we will have to copy it anyway; taking
    // the argument by value forces the copy at the point of the call, and
    // allows passing an rvalue.
    Interp2D(std::vector<Point3D> data);

    // Destructor must clean up allocated GSL resources.
    ~Interp2D() noexcept;

    // Interp2D objects can not be copied because the GSL resources can not
    // be copied.
    Interp2D(Interp2D const&) = delete;

    double operator()(double x, double y) const;

    double
    eval(double x, double y) const
    {
      return this->operator()(x, y);
    };

    // Return the number of grid points in x and y.
    std::size_t nx() const;
    std::size_t ny() const;

  private:
    std::vector<double> xs_;
    std::vector<double> ys_;
    std::vector<double> zs_; // stores z-values in row-major
    gsl_interp2d* interp_;

    // Discover the (x,y) grid implicit in the supplied set of points, or throw
    // a std::domain_error if no grid can be constructed from these points.
    void make_grid_(std::vector<Point3D>& data);
  };
}

template <std::size_t M, std::size_t N>
inline y3_cluster::Interp2D::Interp2D(
  std::array<double, M> const& xs,
  std::array<double, N> const& ys,
  std::array<std::array<double, M>, N> const& zs)
  : xs_(begin(xs), end(xs))
  , ys_(begin(ys), end(ys))
  , zs_(&zs[0][0], &zs[M - 1][N - 1])
{}

inline std::size_t
y3_cluster::Interp2D::nx() const
{
  return xs_.size();
}

inline std::size_t
y3_cluster::Interp2D::ny() const
{
  return ys_.size();
}
#endif
