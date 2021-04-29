#include "catch2/catch.hpp"
#include "utils/make_grid_points.hh"

#include <array>
#include <string>
#include <vector>

std::string const module_label("something");

TEST_CASE("make_grid_splatted")
{

  std::vector<double> radii{1.5, 2.5};
  std::vector<double> zs{3.5};
  std::vector<std::string> names{"radii", "zs"};

  y3_cluster::grid_t<2> mygrid =
    y3_cluster::detail::make_grid_splatted(names, radii, zs);

  CHECK(mygrid.names == names);
  for (std::size_t i = 0; i != mygrid.names.size(); ++i) {
    CHECK(mygrid.points[i][0] == radii[i]); // radius is first coordinate
    CHECK(mygrid.points[i][1] == zs[0]);    // z is second coordinate
  }
}