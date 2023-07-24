#include "cuda_interp_2d.cuh"
#include <utility> // for std::swap

gpu_support::Interp2D::Interp2D(std::vector<double> const& xs,
                                std::vector<double> const& ys,
                                std::vector<double> const& zs) :
  on_host_(xs, ys, zs),
  on_device_(xs, ys, zs)
{ }

gpu_support::Interp2D::Interp2D(std::vector<y3_cluster::Point3D> const& data)
    : on_host_(data)
    , on_device_()
{ }

void
gpu_support::Interp2D::swap(Interp2D& other) noexcept
{
  using std::swap;
  swap(on_host_, other.on_host_);
  swap(on_device_, other.on_device_);
}



