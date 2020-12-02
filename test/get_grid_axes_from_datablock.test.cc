#include "catch2/catch.hpp"
#include "cosmosis/datablock/datablock.hh"
#include "utils/make_grid_points.hh"

#include <array>
#include <string>
#include <vector>

std::string const module_label("something");

TEST_CASE("get_grid_axes_from_datalbock")
{
  cosmosis::DataBlock cfg;
  cfg.put_val(module_label, "radii", std::vector<double>{3.0, 4.0});
  cfg.put_val(module_label, "zs", std::vector<double>{2.5});

  std::array<std::string, 2> names {"radii", "zs"};
  std::array<std::vector<double>, 2> axes;
  y3_cluster::detail::get_grid_axes_from_datablock(
    cfg, module_label, names, axes);
  CHECK(names[0] == "radii");
  CHECK(names[1] == "zs");
}