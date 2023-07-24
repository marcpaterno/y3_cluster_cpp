#ifndef Y3_CLUSTER_CPP_UTILS_CUDA_INTERP_1D_CUH
#define Y3_CLUSTER_CPP_UTILS_CUDA_INTERP_1D_CUH

#include <array>
#include <vector>

#include "interp_1d.hh"
#include "common/cuda/Interp1D.cuh"

// Interp1D is used for linear interpolation in 1 dimension.
// Use this is CUDA models you want to be able to test on a CPU. And if you're
// writing a CUDA model, you should want to test it. And that means testing it
// on a CPU, at least until we know how to make good testing infrastrcture for
// use on the GPU.

namespace gpu_support {

  class Interp1D {
  public:
    // Interpolator created from two arrays; compiler assures they are of the
    // same length.
    template <size_t N>
    Interp1D(std::array<double, N> const& xs, std::array<double, N> const& ys);

    // Interpolator created from two vectors; size equality is tested at
    // runtime.
    Interp1D(std::vector<double> const& xs, std::vector<double> const& ys);

    void swap(Interp1D& other) noexcept;

    __host__ __device__
    double operator()(double x) const;

    __host__ __device__
    double eval(double x) const;

    __host__ __device__
    double clamp(double x) const;

  private:
    // Interp1D contains both a host (CPU) and device (GPU) interpolation
    // table. Both are initialized upon construction, and kept up-to-date
    // in assignment, swapping, etc.
    // When operator()() is called from host code, the host interpolator
    // is used. When operator()() is called from device code, the device
    // interpolator is used.
    y3_cluster::Interp1D on_host_;
    quad::Interp1D on_device_;
  };

} // namespace gpu_support

template <size_t N>
gpu_support::Interp1D::Interp1D(std::array<double, N> const& xs,
                                std::array<double, N> const& ys)
  : on_host_(xs, ys)
  , on_device_(xs, ys)
{ }

__host__ __device__ inline double
gpu_support::Interp1D::operator()(double x) const
{
#ifdef __CUDA_ARCH__
  return on_device_(x);
#else
  return on_host_(x);
#endif
}

__host__ __device__ inline double
gpu_support::Interp1D::eval(double x) const
{
  return (*this)(x);
};

__host__ __device__ inline double
gpu_support::Interp1D::clamp(double x) const
{
#ifdef __CUDA_ARCH__
  return on_device_.clamp(x);
#else
  return on_host_.clamp(x);
#endif
}

#endif
