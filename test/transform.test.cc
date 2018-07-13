#include "catch2/catch.hpp"
#include "transform.hh"
#include <array>
#include <locale> // for std::toupper
#include <string>
#include <type_traits>

TEST_CASE("empty pod")
{
  std::array<int, 0> const xs;
  auto const ys = y3_cluster::transform(xs, [](double x) { return x; });
  CHECK(ys.size() == 0);
}

TEST_CASE("nonempty pod constexpr")
{
  std::array<double, 2> const constexpr xs = {1.5, 2.5};
  auto const constexpr ys =
    y3_cluster::transform(xs, [](double x) { return 2. * x; });
  static_assert(ys[0] == 3.0);
  static_assert(ys[1] == 5.0);
  CHECK(ys[0] == 3.0);
  CHECK(ys[1] == 5.0);
}

TEST_CASE("nonempty pod")
{
  std::array<double, 2> const xs = {1.5, 2.5};
  auto const ys = y3_cluster::transform(xs, [](double x) { return 2. * x; });
  CHECK(ys[0] == 3.0);
  CHECK(ys[1] == 5.0);
}

TEST_CASE("empty nonpod")
{
  std::array<std::string, 0> const xs;
  auto const ys =
    y3_cluster::transform(xs, [](std::string const&) { return std::string(); });
  CHECK(ys.size() == 0);
}

std::string
make_upper(std::string const& s)
{
  std::string result(s);
  std::locale dflt;
  for (auto& ch : result) {
    ch = std::toupper(ch, dflt);
  }
  return result;
}

TEST_CASE("nonempty nonpod")
{
  std::array<std::string, 2> const xs{"cow", "dog"};
  auto const ys = y3_cluster::transform(xs, make_upper);
  CHECK(ys[0] == "COW");
  CHECK(ys[1] == "DOG");
}

struct NDC {
  NDC() = delete;
  int i;
};

TEST_CASE("nonemtpy non-default-constructible")
{
  static_assert(std::is_default_constructible<NDC>::value == false);
  std::array<NDC, 3> xs{-1, 0, 1};
  auto const ys = y3_cluster::transform(xs, [](NDC x) { return x; });
  CHECK(ys[0].i == -1);
  CHECK(ys[1].i == 0);
  CHECK(ys[2].i == 1);
}
