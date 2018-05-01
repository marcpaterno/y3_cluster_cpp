#include "test/interp_2d.hh"
#include "test/fpsupport.hh"

#include <algorithm>
#include <stdexcept>

y3_cluster::Interp2D::Interp2D(std::vector<Point3D> data)
  : xs_(), ys_(), zs_(), interp_(nullptr)
{
  make_grid_(data);
  interp_ = gsl_interp2d_alloc(gsl_interp2d_bilinear, nx(), ny());
  gsl_interp2d_init(interp_, xs_.data(), ys_.data(), zs_.data(), nx(), ny());
}

double
y3_cluster::Interp2D::operator()(double x, double y) const
{
  // We do not use the accelerator features of GSL interpolation, because we
  // do not expect that the pattern of calls will be such that it will help.
  // Profile the resulting integration routine to see if this should be changed.
  return gsl_interp2d_eval_extrap(
    interp_, xs_.data(), ys_.data(), zs_.data(), x, y, nullptr, nullptr);
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
  // This is determining the value of N, the length of a row, which is the
  // number of columns.
  auto equal = [reltol, abstol, pa = data.begin()](Point3D const& b) {
    return fpsupport::is_equivalent((*pa)[0], b[0], reltol, abstol);
  };

  auto const row_1_end = std::find_if_not(data.begin(), data.end(), equal);
  if (row_1_end == data.end())
    throw std::domain_error("Only one row");

  auto const candidate_ncolumns = std::distance(data.begin(), row_1_end);

  // Make sure the number of points is appropriate; determine the number of
  // rows.
  auto const candidate_M = data.size() / candidate_ncolumns;
  if (candidate_M * candidate_ncolumns != data.size())
    throw std::domain_error("Not an integral number of rows");

  // Make sure remaining rows conform: right number of columns, matching
  // y values.
  auto const p_firstrow = data.cbegin();
  auto equivalent_y_values = [reltol, abstol](Point3D const& a,
                                              Point3D const& b) {
    return fpsupport::is_equivalent(a[1], b[1], reltol, abstol);
  };
  for (std::size_t row = 1; row != candidate_M; ++row) {
    auto res = std::mismatch(p_firstrow,
                             p_firstrow + candidate_ncolumns,
                             p_firstrow + row * candidate_ncolumns,
                             equivalent_y_values);
    if (res.first != p_firstrow + candidate_ncolumns)
      throw std::domain_error("Points do not form a grid");
  }
  // Capture x-, y-arrays, and z-matrix.
  throw std::logic_error("This function is not yet complete.");
}
