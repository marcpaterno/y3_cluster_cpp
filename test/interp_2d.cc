#include "test/interp_2d.hh"
#include "test/fpsupport.hh"

#include <algorithm>
#include <stdexcept>

double
y3_cluster::Interp2D::operator()(double x, double y) const
{
  // We do not use the accelerator features of GSL interpolation, because we
  // do not expect that the pattern of calls will be such that it will help.
  // Profile the resulting integration routine to see if this should be changed.
  return gsl_interp2d_eval_extrap(
    interp_, xs_.data(), ys_.data(), zs_.data(), x, y, nullptr, nullptr);
}

y3_cluster::Interp2D::Interp2D(std::vector<Point3D> data)
  : xs_(), ys_(), interp_(nullptr)
{
  make_grid_(data);
  interp_ = gsl_interp2d_alloc(gsl_interp2d_bilinear, nx(), ny());
  gsl_interp2d_init(interp_, xs_.data(), ys_.data(), zs_.data(), nx(), ny());
}

y3_cluster::Interp2D::~Interp2D() noexcept
{
  gsl_interp2d_free(interp_);
}

void
y3_cluster::Interp2D::make_grid_(std::vector<Point3D>& data)
{
  // Discover the (x,y) grid implied by the points we are given -- or reject
  // them as not a grid. Points in x are equivalent if they differ by an
  // absolute value less than xfuzz, and similarly for points in y.
  assure_clean_floats(data);

  // Sort first by x, then y, then z, within equivalence classes determined
  // by given relative and absolute tolerances.
  // TODO: Determine what tolerances are actually useful, and how they should
  // be exposed in the interface.
  double const reltol = 1.e-6;
  double const abstol = 1.e-24;
  Point3DLess comp{reltol, abstol, reltol, abstol};
  std::sort(data.begin(), data.end(), comp);

  // determine length of first block of points (how many y values for that x).
  auto not_equal = [reltol, abstol](Point3D const& a, Point3D const& b){
    return fpsupport::is_equivalent(a[0], b[0], reltol, abstol);
  };
  auto block_1_end = std::adjacent_find(data.begin(), data.end(), not_equal);
  // Make sure remaining blocks all carry appropriately matching y values.
  // Capture x-, y-arrays, and z-matrix.
}
