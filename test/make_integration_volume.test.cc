#include "catch2/catch.hpp"
#include "cosmosis/datablock/datablock.hh"
#include "utils/make_integration_volumes.hh"

#include <string>
#include <vector>

std::string const module_label("something");

using cubacpp::IntegrationVolume;

TEST_CASE("2d volume")
{
  SECTION("vector length 1")
  {
    cosmosis::DataBlock cfg;
    cfg.put_val(module_label, "x_low", std::vector<double>{-3.0});
    cfg.put_val(module_label, "x_high", std::vector<double>{2.0});
    cfg.put_val(module_label, "yyyy_low", std::vector<double>{4.0});
    cfg.put_val(module_label, "yyyy_high", std::vector<double>{7.0});
    auto volumes = y3_cluster::make_integration_volumes_wall_of_numbers(
      cfg, module_label, "x", "yyyy");
    CHECK(volumes.size() == 1);
    IntegrationVolume<2> expected({-3.0, 4.0}, {2.0, 7.0});
    CHECK(volumes[0] == expected);
    CHECK(volumes[0].jacobian() == 5.0 * 3.0);
  }
  SECTION("vector length 2")
  {
    cosmosis::DataBlock cfg;
    cfg.put_val(module_label, "x_low", std::vector<double>{0.0, -3.0});
    cfg.put_val(module_label, "x_high", std::vector<double>{1.0, 2.0});
    cfg.put_val(module_label, "y_low", std::vector<double>{0.0, 4.0});
    cfg.put_val(module_label, "y_high", std::vector<double>{1.0, 7.0});
    auto volumes = y3_cluster::make_integration_volumes_wall_of_numbers(
      cfg, module_label, "x", "y");
    CHECK(volumes.size() == 2);
    IntegrationVolume<2> expected({-3.0, 4.0}, {2.0, 7.0});
    CHECK(volumes[0] == IntegrationVolume<2>{});
    CHECK(volumes[1] == expected);
  }
  SECTION("cartesian product")
  {
    cosmosis::DataBlock cfg;
    cfg.put_val(module_label, "x_low", std::vector<double>{-1.0, 0.0, 1.0});
    cfg.put_val(module_label, "x_high", std::vector<double>{1.0, 2.0, 3.0});
    cfg.put_val(module_label, "y_low", std::vector<double>{0.0, 10.0});
    cfg.put_val(module_label, "y_high", std::vector<double>{10.0, 20.});
    auto volumes = y3_cluster::make_integration_volumes_cartesian_product(
      cfg, module_label, "x", "y");
    CHECK(volumes.size() == 6);
    std::vector<IntegrationVolume<2>> expected;
    expected.reserve(6);
    expected.push_back(IntegrationVolume<2>({-1.0, 0.0}, {1.0, 10.0}));
    expected.push_back(IntegrationVolume<2>({-1.0, 10.0}, {1.0, 20.0}));
    expected.push_back(IntegrationVolume<2>({0.0, 0.0}, {2.0, 10.0}));
    expected.push_back(IntegrationVolume<2>({0.0, 10.0}, {2.0, 20.0}));
    expected.push_back(IntegrationVolume<2>({1.0, 0.0}, {3.0, 10.0}));
    expected.push_back(IntegrationVolume<2>({1.0, 10.0}, {3.0, 20.0}));
    for (std::size_t i = 0; i != 6; ++i) { CHECK(volumes[i] == expected[i]); }
  }
}
