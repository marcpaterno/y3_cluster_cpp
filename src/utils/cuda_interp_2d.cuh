#ifndef Y3_CLUSTER_CPP_UTILS_CUDA_INTERP_2D_CUH
#define Y3_CLUSTER_CPP_UTILS_CUDA_INTERP_2D_CUH

#include <array>
#include <vector>

#include "interp_2d.hh"
#include "point_3d.hh"
#include "common/cuda/Interp2D.cuh"

namespace gpu_support {

  class Interp2D {
    public:
      // We publish Interp2D::matrix as part of our public interface because some of
      // the testing on CPU uses this.
      template <size_t M, size_t N> using matrix = y3_cluster::Interp2D::matrix<M, N>;

      // Interpolator created from 3 arrays; compiler assures that zs matches
      // the dimensionality of xs and ys.
      template <size_t M, size_t N>
      Interp2D(std::array<double, M> const& xs,
               std::array<double, N> const& ys,
               matrix<M,N> const& zs);

      // Interpolator create from three vectors; sizes are checked for
      // consistency at runtime.
      Interp2D(std::vector<double> const& xs,
               std::vector<double> const& ys,
               std::vector<double> const& zs);

      // Interpolator created from two vectors and one vector-of-vectors.
      Interp2D(std::vector<double> const& xs,
               std::vector<double> const& ys,
               std::vector<std::vector<double>> const& zs);

      // Like above - assumes ndarray is 2D and fits it appropriately
      Interp2D(std::vector<double> const& xs,
               std::vector<double> const& ys,
               cosmosis::ndarray<double> const& zs);

      // Create only the CPU-based interpolator (and not the GPU-based
      // interpolator) from vector of triplets. Throws std::logic_error if
      // the points do not lie on a grid in (x,y) space, or if any values are
      // NaN or infinities. Any denormalized x- or y-values are flushed to
      // zero.
      //
      // It is an error to attempt to use this object on the GPU; such use is
      // not diagnosed at runtime, and will probably result in a difficult-to-
      // diagnose crash.
      Interp2D(std::vector<y3_cluster::Point3D> const& data);

      void swap(Interp2D& other) noexcept;

      __host__ __device__
      double operator()(double x, double y) const;

      __host__ __device__
      double eval(double x, double y) const;

      __host__ __device__
      double clamp(double x, double y) const;

    private:

    // Interp2D contains both a host (CPU) and device (GPU) interpolation
    // table. Both are initialized upon construction, and kept up-to-date
    // in assignment, swapping, etc.
    // When operator()() is called from host code, the host interpolator
    // is used. When operator()() is called from device code, the device
    // interpolator is used.
    y3_cluster::Interp2D on_host_;
    quad::Interp2D on_device_;
  };

} // namespace gpu_support

template <size_t M, size_t N>
gpu_support::Interp2D::Interp2D(std::array<double, M> const& xs,
                                std::array<double, N> const& ys,
                                matrix<M, N> const& zs)
  : on_host_(xs, ys, zs)
  , on_device_(xs, ys, zs)
{ }

__host__ __device__ inline double
gpu_support::Interp2D::operator()(double x, double y) const
{
#ifdef __CUDA_ARCH__
  return on_device_(x, y);
#else
  return on_host_(x, y);
#endif
}

__host__ __device__ inline double
gpu_support::Interp2D::eval(double x, double y) const
{
  return (*this)(x, y);
}

__host__ __device__ inline double
gpu_support::Interp2D::clamp(double x, double y) const
{
#ifdef __CUDA_ARCH__
  return on_device_.clamp(x, y);
#else
  return on_host_.clamp(x, y);
#endif
}

#endif
