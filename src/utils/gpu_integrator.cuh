#ifndef Y3_CLUSTER_UTIL_CUDA_INTEGRATOR_CUH
#define Y3_CLUSTER_UTIL_CUDA_INTEGRATOR_CUH

#include "cubacpp/arity.hh"
#include "cuda/cudaPagani/quad/GPUquad/Pagani.cuh"
#include "cuda/cudaPagani/quad/util/Volume.cuh"
#include "cuda/mcubes/mcubes.cuh"
#include "cuda/mcubes/vegasT.cuh"
#include "cuda/mcubes/vegasSeqMcubes.hh"

#include <tuple>

namespace y3_cuda {
  template <int ndim>
  class MultiDimensionalIntegrator {
  public:
    std::size_t num_algs() const;

    MultiDimensionalIntegrator() = default;
    explicit MultiDimensionalIntegrator(std::string const& name);

    template <class F>
    cuhreResult<double> integrate(int which,
                                  F f,
                                  double epsabs,
                                  double epsrel);

    template <class F>
    cuhreResult<double> integrate(int which,
                                  F f,
                                  double epsabs,
                                  double epsrel,
                                  quad::Volume<double, ndim> const* vol);

    template <class F>
    cuhreResult<double> integrate(F f, double epsabs, double epsrel);

    template <class F>
    cuhreResult<double> integrate(F f,
                                  double epsabs,
                                  double epsrel,
                                  quad::Volume<double, ndim> const* vol);

    void set_maxeval(double m);

  private:
    using algs_t = std::tuple<quad::Pagani<double, ndim>, quad::mcubes, VegasSEQmcubes>;
    algs_t algorithms_;
    int which_ = 0;
  };
}

template <int N>
y3_cuda::MultiDimensionalIntegrator<N>::MultiDimensionalIntegrator(
  std::string const& name)
{
  if (name == std::string("pagani"))
    which_ = 0;
  else if (name == std::string("mcubes"))
    which_ = 1;
  else if(name == std::string("seqmcubes"))
    which_ = 2;
  else
    throw std::runtime_error("MultiDimensionalIntegrator::integrate called for "
                             "unrecognized algorithm");
}

template <int N>
std::size_t
y3_cuda::MultiDimensionalIntegrator<N>::num_algs() const
{
  return std::tuple_size_v<algs_t>;
}

template <int N>
template <class F>
cuhreResult<double>
y3_cuda::MultiDimensionalIntegrator<N>::integrate(int which,
                                               F f,
                                               double epsabs,
                                               double epsrel)
{
  switch (which) {
    case 0:
      return std::get<0>(algorithms_)
        .integrate(std::forward<F>(f), epsabs, epsrel);
    case 1:
      return std::get<1>(algorithms_)
        .integrate(std::forward<F>(f), epsabs, epsrel);
    case 2:
      return std::get<3>(algorithms_)
        .integrate(std::forward<F>(f), epsabs, epsrel);
  default:
      throw std::runtime_error("MultiDimensionalIntegrator::integrate called "
                               "for unrecognized algorithm");
  }
}

template <int N>
template <class F>
cuhreResult<double>
y3_cuda::MultiDimensionalIntegrator<N>::integrate(
  int which,
  F f,
  double epsabs,
  double epsrel,
  quad::Volume<double, N> const* vol)
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

template <int N>
template <class F>
cuhreResult<double>
y3_cuda::MultiDimensionalIntegrator<N>::integrate(F f,
                                               double epsabs,
                                               double epsrel)
{
  switch (which_) {
    case 0:
      return std::get<0>(algorithms_)
        .integrate(std::forward<F>(f), epsabs, epsrel);
    case 1:
      return std::get<1>(algorithms_)
        .integrate(std::forward<F>(f), epsabs, epsrel);
    case 2:
      return std::get<3>(algorithms_)
        .integrate(std::forward<F>(f), epsabs, epsrel);
    
    default:
      throw std::runtime_error("MultiDimensionalIntegrator::integrate called "
                               "for unrecognized algorithm");
  }
}

template <int N>
template <class F>
cuhreResult<double>
y3_cuda::MultiDimensionalIntegrator<N>::integrate(
  F f,
  double epsabs,
  double epsrel,
  quad::Volume<double, N> const* vol)
{
  switch (which_) {
    case 0:
      return std::get<0>(algorithms_)
        .integrate(f, epsrel, epsabs, vol);
  case 1:{
      auto x = std::get<1>(algorithms_)
        .integrate(f, epsabs, epsrel, vol);
      return x;
      break;
  }

  case 2:{
    auto x = std::get<2>(algorithms_)
        .integrate(f, epsabs, epsrel, vol);
      return x;
      break;
  }
    default:
      throw std::runtime_error("MultiDimensionalIntegrator::integrate called "
                               "for unrecognized algorithm");
  }
}

template <int N>
void
y3_cuda::MultiDimensionalIntegrator<N>::set_maxeval(double m)
{
  if (m <= 0.0) {
    throw std::runtime_error("MultiDimensionalIntegrator: maxcalls can not be negative\n");
  }

  // PAGANI does not use maxcalls
  std::get<1>(algorithms_).maxcalls = static_cast<long long>(m);
  std::get<2>(algorithms_).maxcalls = static_cast<unsigned long>(m);
}

#endif
