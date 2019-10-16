#ifndef __Y3_PARAM_SPACE_EXPLORER
#define __Y3_PARAM_SPACE_EXPLORER

#include <algorithm>
#include <array>
#include <cstddef>
#include <numeric>

#include "utils/integration_range.hh"

/* Explore a grid of parameters.
 */
template <std::size_t NParams>
struct ParamSpaceExplorer {
  template <typename T>
  using array = std::array<T, NParams>;
  array<double*> values;
  array<const y3_cluster::IntegrationRange> ranges;
  array<const std::size_t> width;
  array<std::size_t> cur_pos;
  std::size_t iterations, max_iterations;

  ParamSpaceExplorer(array<double*> values_,
                     array<const y3_cluster::IntegrationRange> ranges_,
                     array<const std::size_t> width_)
    : values(values_), ranges(ranges_), width(width_), iterations(0)
  {
    max_iterations = 1;
    for (auto i = 0u; i < NParams; i++) {
      cur_pos[i] = 0;
      *values[i] = ranges_[i].transform(0.0);
      max_iterations *= width[i];
    }
  }

  bool
  next()
  {
    if (++iterations >= max_iterations)
      return false;

    std::size_t i{NParams - 1};

    while (cur_pos[i] == (width[i] - 1)) {
      *values[i] = ranges[i].transform(0.0);
      cur_pos[i--] = 0;
    }

    *values[i] = ranges[i].transform(++cur_pos[i] / (width[i] - 1.0));

    return true;
  }
};

#endif
