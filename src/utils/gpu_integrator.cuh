#ifndef Y3_CLUSTER_UTIL_CUDA_INTEGRATOR_CUH
#define Y3_CLUSTER_UTIL_CUDA_INTEGRATOR_CUH

#include "cudaPagani/quad/GPUquad/Pagani.cuh"
#include "vegas/vegasT.cuh"

#include <tuple>

namespace y3_cuda {
  class MultiDimensionalIntegrator {
  public:
    std::size_t num_algs() const;

    MultiDimensionalIntegrator() = default;
    explicit MultiDimensionalIntegrator(std::string const& name);

    template <class F>
    cubacpp::integration_result integrate(int which,
                                          F f,
                                          double epsrel,
                                          double epsabs) const;

    template <class F>
    cubacpp::integration_result integrate(
      int which,
      F f,
      double epsrel,
      double epsabs,
      typename cubacpp::integration_volume_for<F>::type vol) const;

    template <class F>
    cubacpp::integration_result integrate(F f,
                                          double epsrel,
                                          double epsabs) const;

    template <class F>
    cubacpp::integration_result integrate(
      F f,
      double epsrel,
      double epsabs,
      typename cubacpp::integration_volume_for<F>::type vol) const;

    void set_maxeval(double m);

  private:
    using algs_t = std::tuple<quad::Pagani, quad::mcubes>;
    algs_t algorithms_;
    int which_ = 0;
  };
}

inline y3_cuda::MultiDimensionalIntegrator::MultiDimensionalIntegrator(
  std::string const& name)
{
  if (name == std::string("pagani"))
    which_ = 0;
  else if (name == std::string("mcubes"))
    which_ = 1;
  else
    throw std::runtime_error("MultiDimensionalIntegrator::integrate called for "
                             "unrecognized algorithm");
}

inline std::size_t
y3_cuda::MultiDimensionalIntegrator::num_algs() const
{
  return std::tuple_size_v<algs_t>;
}

template <class F>
cubacpp::integration_result
y3_cuda::MultiDimensionalIntegrator::integrate(int which,
                                                  F f,
                                                  double epsabs,
                                                  double epsrel) const
{
  switch (which) {
    case 0:
      return std::get<0>(algorithms_)
        .integrate(std::forward<F>(f), epsabs, epsrel);
    case 1:
      return std::get<1>(algorithms_)
        .integrate(std::forward<F>(f), epsabs, epsrel);
    case 2:
      return std::get<2>(algorithms_)
        .integrate(std::forward<F>(f), epsabs, epsrel);
    default:
      throw std::runtime_error("MultiDimensionalIntegrator::integrate called "
                               "for unrecognized algorithm");
  }
}

template <class F>
cubacpp::integration_result
y3_cuda::MultiDimensionalIntegrator::integrate(
  int which,
  F f,
  double epsabs,
  double epsrel,
  typename cubacpp::integration_volume_for<F>::type vol) const
{
  switch (which) {
    case 0:
      return std::get<0>(algorithms_)
        .integrate(std::forward<F>(f), epsabs, epsrel, vol);
    case 1:
      return std::get<1>(algorithms_)
        .integrate(std::forward<F>(f), epsabs, epsrel, vol);
    case 2:
      return std::get<2>(algorithms_)
        .integrate(std::forward<F>(f), epsabs, epsrel, vol);
    default:
      throw std::runtime_error("MultiDimensionalIntegrator::integrate called "
                               "for unrecognized algorithm");
  }
}

template <class F>
cubacpp::integration_result
y3_cuda::MultiDimensionalIntegrator::integrate(F f,
                                                  double epsabs,
                                                  double epsrel) const
{
  switch (which_) {
    case 0:
      return std::get<0>(algorithms_)
        .integrate(std::forward<F>(f), epsabs, epsrel);
    case 1:
      return std::get<1>(algorithms_)
        .integrate(std::forward<F>(f), epsabs, epsrel);
    case 2:
      return std::get<2>(algorithms_)
        .integrate(std::forward<F>(f), epsabs, epsrel);
    default:
      throw std::runtime_error("MultiDimensionalIntegrator::integrate called "
                               "for unrecognized algorithm");
  }
}

template <class F>
cubacpp::integration_result
y3_cuda::MultiDimensionalIntegrator::integrate(
  F f,
  double epsabs,
  double epsrel,
  typename cubacpp::integration_volume_for<F>::type vol) const
{
  switch (which_) {
    case 0:
      return std::get<0>(algorithms_)
        .integrate(std::forward<F>(f), epsabs, epsrel, vol);
    case 1:
      return std::get<1>(algorithms_)
        .integrate(std::forward<F>(f), epsabs, epsrel, vol);
    case 2:
      return std::get<2>(algorithms_)
        .integrate(std::forward<F>(f), epsabs, epsrel, vol);
    default:
      throw std::runtime_error("MultiDimensionalIntegrator::integrate called "
                               "for unrecognized algorithm");
  }
}

void
y3_cuda::MultiDimensionalIntegrator::set_maxeval(long long int m)
{
  std::get<0>(algorithms_).maxeval = m;
  std::get<1>(algorithms_).maxeval = m;
}

#endif
