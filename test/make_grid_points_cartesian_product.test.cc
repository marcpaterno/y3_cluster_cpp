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

  CHECK(mygrid.names.size() == 2UL);
  CHECK(mygrid.names[0] == "zs");
  CHECK(mygrid.names[1] == "radii");
  CHECK(mygrid.points.size() == 2UL);
  CHECK(mygrid.points[0] == std::array<double, 2>{2.5, 3.0});
  CHECK(mygrid.points[1] == std::array<double, 2>{2.5, 4.0});
}

TEST_CASE("2 x 5 x 3 grid")
{
  cosmosis::DataBlock cfg;
  cfg.put_val(module_label, "zo_low", std::vector<double>{-1.5, 1.5});
  cfg.put_val(module_label, "zo_high", std::vector<double>{1.25, 2.25, 3.25, 4.25, 5.25});
  cfg.put_val(module_label, "rmis_array", std::vector<double>{1., 2., 4.});
  y3_cluster::grid_t mygrid = y3_cluster::make_grid_points_cartesian_product(cfg, module_label, "zo_low", "zo_high", "rmis_array");
  CHECK(mygrid.names.size() == 3UL);
  CHECK(mygrid.names[0] == "zo_low");
  CHECK(mygrid.names[1] == "zo_high");
  CHECK(mygrid.names[2] == "rmis_array");
  CHECK(mygrid.points.size() == 30UL);
  CHECK(mygrid.points[0] == std::array<double, 3>{-1.5, 1.25, 1.});
  CHECK(mygrid.points[29] == std::array<double, 3>{1.5, 5.25, 4.});
}
