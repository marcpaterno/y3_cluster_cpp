#ifndef Y3_CLUSTER_CPP_MAKE_INTEGRATION_VOLUME_HH
#define Y3_CLUSTER_CPP_MAKE_INTEGRATION_VOLUME_HH

#include "cosmosis/datablock/datablock.hh"
#include "cubacpp/array.hh"
#include "cubacpp/integration_volume.hh"

#include <array>
#include <string>
#include <type_traits>
#include <vector>

namespace y3_cluster {
  // This variadic function template takes as arguments:
  //   1. a cosmosis::DataBlock (by reference), and
  //   2. one or more arguments that are convertible to strings.
  //
  // It should be used like:
  //    DataBlock cfg;
  //    auto vols = make_integration_volumes(cfg, "x", "y", "z");
  //
  // The function returns a vector of IntegrationVolume<N>,
  // where N is the number of names provided.
  template <typename... Ts>
  std::vector<cubacpp::IntegrationVolume<sizeof...(Ts)>>
  make_integration_volumes(cosmosis::DataBlock& cfg,
                           std::string const& modulelabel,
                           Ts... names);

  // We have to use int N, rather than std::size_t N, because that is what the
  // class template cubacpp::array<N> expects; this in turn is determined by
  // Eigen.
  template <std::size_t N>
  void
  get_integration_boundaries(cosmosis::DataBlock& cfg,
                             std::string const& modulelabel,
                             std::array<std::string, N> const& names,
                             std::vector<cubacpp::array<N>>& lowbounds,
                             std::vector<cubacpp::array<N>>& highbounds)
  {
    if (names.empty())
      return;

    // The first parameter gets special handling, because we use it to determine
    // how many integration volumes we shall produce.
    using vec = std::vector<double>;
    auto lows = cfg.view<vec>(modulelabel, names[0] + "_low");
    std::size_t const nvolumes = lows.size();

    auto highs = cfg.view<vec>(modulelabel, names[0] + "_high");
    if (nvolumes != highs.size()) {
      // TODO: Improve this error handling.
      throw std::runtime_error(names[0]+" bad, bad user!");
    }

    lowbounds.resize(nvolumes);
    highbounds.resize(nvolumes);
    auto fill_bounds =
      [&lowbounds, &lows, &highbounds, &highs](std::size_t iname, int ivol) {
        lowbounds[ivol][iname] = lows[ivol];
        highbounds[ivol][iname] = highs[ivol];
      };

    for (std::size_t ivol = 0; ivol != nvolumes; ++ivol)
      fill_bounds(0, ivol);

    // All other parameters are handled identically to each other. Each must
    // have the same number of integration volumes.
    for (std::size_t iname = 1; iname != N; ++iname) {
      lows = cfg.view<vec>(modulelabel, names[iname] + "_low");
      highs = cfg.view<vec>(modulelabel, names[iname] + "_high");
      if (nvolumes != lows.size()) {
        // TODO: Improve this error handling.
        throw std::runtime_error(names[iname]+" bad, bad user!");
      }

      if (nvolumes != highs.size()) {
        abort();
        throw std::runtime_error(names[iname]+" bad, bad user!");
      }
      for (std::size_t ivol = 0; ivol != nvolumes; ++ivol)
        fill_bounds(iname, ivol);
    }
  }
} // namespace y3_cluster

template <typename... Ts>
std::vector<cubacpp::IntegrationVolume<sizeof...(Ts)>>
y3_cluster::make_integration_volumes(cosmosis::DataBlock& cfg,
                                         std::string const& modulelabel,
                                         Ts... ts)
{
  // Make sure that all arguments are convertible to std::string.
  static_assert(std::conjunction_v<std::is_convertible<Ts, std::string>...>,
                "\n\nCosmoSIS error!\nAll trailing arguments in "
                "make_integration_volumes must be convertible to string.\n\n");

  constexpr std::size_t n = sizeof...(Ts);
  std::array<std::string, n> names{std::forward<Ts>(ts)...};

  using boundary_t = cubacpp::array<n>;
  std::vector<boundary_t> lows, highs;
  y3_cluster::get_integration_boundaries(
    cfg, modulelabel, names, lows, highs);

  std::vector<cubacpp::IntegrationVolume<n>> result;
  result.reserve(lows.size());
  for (std::size_t i = 0; i != lows.size(); ++i) {
    result.emplace_back(lows[i], highs[i]);
  }

  return result;
}

#endif
