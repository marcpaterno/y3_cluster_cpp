#include "cuda_interp_1d.cuh"
#include <utility> // for std::swap

gpu_support::Interp1D::Interp1D(std::vector<double> const& xs,
                                std::vector<double> const& ys) :
  on_host_(xs, ys),
  on_device_(xs.data(), ys.data(), xs.size())
{ }

void
gpu_support::Interp1D::swap(Interp1D& other) noexcept
{
  using std::swap;
  swap(on_host_, other.on_host_);
  swap(on_device_, other.on_device_);
}



