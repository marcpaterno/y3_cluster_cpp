#include "catch2/catch.hpp"
#include "utils/make_grid_points.hh"

#include <array>
#include <string>
#include <vector>

std::string const module_label("something");
using namespace y3_cluster;

TEST_CASE("make_grid_cartesian_product")
{
  constexpr std::size_t n_dims = 2;
  std::array<std::vector<double>, n_dims> axes = {{{1.5, 2.5}, {3.5}}};
  std::array<std::string, n_dims> names{{"radii", "zs"}};

  grid_t<n_dims> my_grid = detail::make_grid_cartesian_product(axes, names);
  for (std::size_t i = 0; i != n_dims; ++i) {
    CHECK(my_grid.names[i] == names[i]);
    CHECK(my_grid.points[i][0] == axes[0][i]);
    CHECK(my_grid.points[i][1] == axes[1][0]);
  }
}