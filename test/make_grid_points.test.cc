#include "catch2/catch.hpp"
#include "cosmosis/datablock/datablock.hh"
#include "utils/make_grid_points.hh"

#include <array>
#include <cstdlib>
#include <string>
#include <vector>


std::string const module_label("something");

TEST_CASE("1d grid")
{
  using grid_t = y3_cluster::grid_t<1>;

  SECTION("1x1 grid")
  {
    cosmosis::DataBlock cfg;
    cfg.put_val(module_label, "radii", std::vector<double>{3.0});
    auto grid = y3_cluster::make_grid_points_cartesian_product(
      cfg, module_label, "radii");
    CHECK(grid.size() == 1);
    grid_t expected({{{3.0}}},
    std::vector<std::string>{"radii"});
    CHECK(grid == expected);
  }

  SECTION("1x2 grid")
  {
    cosmosis::DataBlock cfg;
    cfg.put_val(module_label, "radii", std::vector<double>{5.0, 10.0});
    auto grid = y3_cluster::make_grid_points_cartesian_product(
      cfg, module_label, "radii");
    CHECK(grid.size() == 2);
    std::vector<std::array<double, 1>> expected({{{5.0}}, {{10.0}}});
    CHECK(grid.points == expected);
  }
}

TEST_CASE("2d grid")
{
  using grid_t = y3_cluster::grid_t<2>;
 
  SECTION("2x1 grid")
  {
    cosmosis::DataBlock cfg;
    cfg.put_val(module_label, "radii", std::vector<double>{3.0, 4.0});
    cfg.put_val(module_label, "zs", std::vector<double>{2.5});
    auto grid = y3_cluster::make_grid_points_cartesian_product(
      cfg, module_label, "radii", "zs");
    CHECK(grid.size() == 2);
    std::vector<std::array<double,2>> expected({{{3.0, 2.5}}, {{4.0, 2.5}}});
    CHECK(grid.points == expected);
  }

  SECTION("2x1 grid wall of numbers")
  {
    cosmosis::DataBlock cfg;
    cfg.put_val(module_label, "radii", std::vector<double>{3.0, 4.0, 5.0});
    cfg.put_val(module_label, "zs", std::vector<double>{1.5, 1.5, 2.5});
    grid_t grid = y3_cluster::make_grid_points_wall_of_numbers(
      cfg, module_label, "radii", "zs");
    CHECK(grid.size() == 3);
    std::vector<std::array<double,2>> expected{{{3.0, 1.5}}, {{4.0, 1.5}}, {{5.0, 2.5}}};
    CHECK(grid.points == expected);
  }

  SECTION("2x1 grid load from file")
  {
    cosmosis::DataBlock cfg;
    std::string const fname =
      std::string(std::getenv("Y3_CLUSTER_CPP_DIR")) +
      "/data/" + "test_grid_file.txt";
    cfg.put_val(module_label, "grid_file", fname);
    grid_t grid = y3_cluster::load_grid_from_file_wall_of_numbers(
      cfg, module_label, "radii", "zs");
    CHECK(grid.size() == 3);
  }
}

TEST_CASE("3d grid")
{
  using grid_t = y3_cluster::grid_t<3>;
  using grid_point_t = grid_t::value_type;

  std::vector<double> xs{1.0, 2.0, 3.0};
  std::vector<double> ys{4.0, 5.0};
  std::vector<double> zs{6.0, 7.0, 8.0, 9.0};

  grid_t grid = y3_cluster::detail::make_grid_splatted(names, xs, ys, zs);
  CHECK(grid.size() == 3 * 2 * 4);

  std::vector<grid_point_t> expected{
    {1, 4, 6}, {1, 4, 7}, {1, 4, 8}, {1, 4, 9}, {1, 5, 6}, {1, 5, 7},
    {1, 5, 8}, {1, 5, 9}, {2, 4, 6}, {2, 4, 7}, {2, 4, 8}, {2, 4, 9},
    {2, 5, 6}, {2, 5, 7}, {2, 5, 8}, {2, 5, 9}, {3, 4, 6}, {3, 4, 7},
    {3, 4, 8}, {3, 4, 9}, {3, 5, 6}, {3, 5, 7}, {3, 5, 8}, {3, 5, 9}};
  CHECK(grid.points == expected);
  CHECK(grid.names == names);
}
