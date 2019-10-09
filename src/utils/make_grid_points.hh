#ifndef Y3_CLUSTER_CPP_MAKE_GRID_POINTS_HH
#define Y3_CLUSTER_CPP_MAKE_GRID_POINTS_HH

#include "cosmosis/datablock/datablock.hh"
#include "utils/meta.hh"

#include <algorithm>
#include <array>
#include <string>
#include <tuple> // for std::apply
#include <vector>

namespace y3_cluster {
  // This variadic function template takes as arguments:
  //   1. a cosmosis::DataBlock (by reference), and
  //   2. one or more arguments that are convertible to strings.
  //
  // It should be used like:
  //    DataBlock cfg;
  //    auto vols = make_grid_points(cfg, "radius", "z");
  //
  // The function returns a vector of std::array<double, N>,
  // where N is the number of names provided.
  template <typename... Ts>
  std::vector<std::array<double, sizeof...(Ts)>> make_grid_points(
    cosmosis::DataBlock& cfg,
    std::string const& modulelabel,
    Ts... names);

  namespace detail {

    // Given a vector of N names of configuration parameters, look up
    // the N vectors of grid locations for each of the Nm axes, and
    // put them into 'axes'.
    template <std::size_t N>
    void
    get_grid_axes(cosmosis::DataBlock& cfg,
                  std::string const& modulelabel,
                  std::array<std::string, N> const& names,
                  std::array<std::vector<double>, N>& axes)
    {
      auto get_axis = [&cfg, &modulelabel](std::string const& name) {
        return cfg.view<std::vector<double>>(modulelabel, name);
      };
      std::transform(names.begin(), names.end(), axes.begin(), get_axis);
    }

  } // namespace detail

  // make_grid_splatted takes N vectors of floating point numbers, where N is
  // any positive integer, and returns a vector of grid points, where each grid
  // point is an array of N doubles.
  template <typename... Ts>
  std::vector<std::array<double, sizeof...(Ts)>>
  make_grid_splatted(std::vector<Ts> const&... axes)
  {
    // Make sure all the vectors are carrying floating point numbers;
    // we convert everything to double, 'cause doing otherwise is hard.
    static_assert(std::conjunction_v<std::is_floating_point<Ts>...>);

    std::vector<std::array<double, sizeof...(Ts)>> res;
    auto accumulator = [&res](Ts const&... ts) {
      // Construct an array from all the elements we pass in ts...
      res.push_back({ts...});
    };
    detail::cartesian_product_aux(accumulator, axes...);
    return res;
  }

  // make_grid is called with an array of arbitrary length, and an
  // index_sequence of matching length.
  template <std::size_t... Is>
  std::vector<std::array<double, sizeof...(Is)>>
  make_grid_aux(std::array<std::vector<double>, sizeof...(Is)> const& axes,
                std::index_sequence<Is...> /*unused*/)
  {
    return make_grid_splatted(axes[Is]...);
  }

  // use this one!
  template <std::size_t N>
  std::vector<std::array<double, N>>
  make_grid(std::array<std::vector<double>, N> const& axes)
  {
    return make_grid_aux(axes, std::make_index_sequence<N>());
  }

} // namespace y3_cluster

// This is the function users of CosmoSIS should use to create a
// grid for the calculation of integrals.
template <typename... STRINGLIKEs>
std::vector<std::array<double, sizeof...(STRINGLIKEs)>>
y3_cluster::make_grid_points(cosmosis::DataBlock& cfg,
                             std::string const& modulelabel,
                             STRINGLIKEs... stringlikes)
{
  // Make sure that all arguments are convertible to std::string.
  static_assert(
    std::conjunction_v<std::is_convertible<STRINGLIKEs, std::string>...>,
    "\n\nCosmoSIS error!\nAll trailing arguments in "
    "make_grid_points must be convertible to string.\n\n");
  constexpr std::size_t n_axes = sizeof...(STRINGLIKEs);

  std::array<std::string, n_axes> const axis_names{
    std::forward<STRINGLIKEs>(stringlikes)...};
  std::array<std::vector<double>, n_axes> axes;
  detail::get_grid_axes(cfg, modulelabel, axis_names, axes);
  return make_grid(axes);
}

#endif
