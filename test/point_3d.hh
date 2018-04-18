#ifndef Y3_CLUSTER_POINT_3D_HH
#define Y3_CLUSTER_POINT_3D_HH

#include <array>

namespace y3_cluster {

  inline bool
  is_equivalent(double x, double y, double relTol, double absTol)
  {
    return (std::abs(x - y) <=
            std::max(absTol, relTol * std::max(std::abs(x), std::abs(y))));
  }

  // Point3D represents an (x,y,z) triplet.
  using Point3D = std::array<double, 3>;

  // Point3DLess is a predicate object for ordering Point3D objects, in a manner
  // suitable for the establishment of (x,y) grids for interpolation in 2D. For
  // identifying the x- and y-coordinates of the grid, approximate equality
  // testing is used. 'xabs' and 'xrel' are the absolute and relative tolerances
  // for evaluating the equality of the x-coordinates, and 'yabs' and 'yrel'
  // those for the y-coordinates.
  class Point3DLess {
  public:
    Point3DLess(double xrel, double xabs, double yrel, double yabs) noexcept;
    bool operator()(Point3D const& a, Point3D const& b) const;

  private:
    double xrel_;
    double xabs_;
    double yrel_;
    double yabs_;
  };
}

inline y3_cluster::Point3DLess::Point3DLess(double xrel,
                                            double xabs,
                                            double yrel,
                                            double yabs) noexcept
  : xrel_(xrel), xabs_(xabs), yrel_(yrel), yabs_(yabs)
{}

inline bool
y3_cluster::Point3DLess::operator()(Point3D const& a, Point3D const& b) const
{
  if (not is_equivalent(a[0], b[0], xabs_, xrel_))
    return a[0] < b[0];
  if (not is_equivalent(a[1], b[1], yabs_, yrel_))
    return a[1] < b[1];
  return a[2] < b[2];
}

#endif