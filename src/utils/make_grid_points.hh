#ifndef Y3_CLUSTER_CPP_MAKE_GRID_POINTS_HH
#define Y3_CLUSTER_CPP_MAKE_GRID_POINTS_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"

#include <algorithm>
#include <array>
#include <string>
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

} // namespace y3_cluster

template <typename... Ts>
std::vector<std::array<double, sizeof...(Ts)>>
y3_cluster::make_grid_points(cosmosis::DataBlock& cfg,
                            std::string const& modulelabel,
                            Ts... ts)
{
  // Make sure that all arguments are convertible to std::string.
  static_assert(std::conjunction_v<std::is_convertible<Ts, std::string>...>,
                "\n\nCosmoSIS error!\nAll trailing arguments in "
                "make_grid_points must be convertible to string.\n\n");
  constexpr std::size_t n_axes = sizeof...(Ts);
  using grid_point = std::array<double, n_axes>;

  std::array<std::string, n_axes> const names{std::forward<Ts>(ts)...};
  std::array<std::vector<double>, n_axes> axes;
  get_grid_axes(cfg, modulelabel, names, axes);
  std::size_t n_grid_points = 1;
  for (auto const& axis : axes) n_grid_points *= axis.size();
  std::vector<grid_point> result(n_grid_points);
  
  for (std::size_t i = 0; i !=  n_axes; ++i) {
    std::size_t const n_points_this_axis = axes[i].size();
    for (std::size_t j = 0; j != n_points_this_axis; ++j) {
      double const current_coor_value = axes[i][j];
      result[i*n_axes + j][i] = current_coor_value;
    }
  }
  return result;
}

#endif
