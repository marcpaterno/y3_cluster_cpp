#include "catch2/catch.hpp"
#include "cosmosis/datablock/datablock.hh"
#include "utils/make_grid_points.hh"

#include <array>
#include <string>
#include <vector>

std::string const module_label("something");

TEST_CASE("make_grid_points_cartesian_product")
{
  cosmosis::DataBlock cfg;
  cfg.put_val(module_label, "radii", std::vector<double>{3.0, 4.0});
  cfg.put_val(module_label, "zs", std::vector<double>{2.5});

  y3_cluster::grid_t mygrid = y3_cluster::make_grid_points_cartesian_product(
    cfg, module_label, "zs", "radii");
  
  CHECK(mygrid.names[0] == "zs");
  CHECK(mygrid.names[1] == "radii");
  CHECK(mygrid.points[0] == std::array<double, 2>{2.5, 3.0});
  CHECK(mygrid.points[1] == std::array<double, 2>{2.5, 4.0});
}