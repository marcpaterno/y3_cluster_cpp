#ifndef Y3_CLUSTER_FPSUPPORT_HH
#define Y3_CLUSTER_FPSUPPORT_HH

#include "test/point_3d.hh"
#include <cmath>

namespace fpsupport {
  // Detect "icky" values: NaN and infinities
  bool
  icky(double x)
  {
    auto const code = std::fpclassify(x);
    return (code == FP_NAN) || (code == FP_INFINITE);
  }

  bool
  icky(y3_cluster::Point3D const& p)
  {
    return icky(p[0]) || icky(p[1]) || icky(p[2]);
  }

  // Replace any subnormal x- or y-values by 0. We assume there are no NaN or
  // infinite values in p.
  void
  squash_subnormals(y3_cluster::Point3D& p)
  {
    if (not std::isnormal(p[0]))
      p[0] = 0.0;
    if (not std::isnormal(p[1]))
      p[1] = 0.0;
  }

  // "Clean" the input points. If any NaN or infinity values are detected,
  // throw std::domain_error. Replace any x- or y-value denormals by zero.
  void
  assure_clean_floats(std::vector<y3_cluster::Point3D>& points)
  {
    for (auto& p : points) {
      if (icky(p))
        throw std::domain_error("Inf. or NaN detected in Interp2D setup");
      squash_subnormals(p);
    }
  }
}

#endif
