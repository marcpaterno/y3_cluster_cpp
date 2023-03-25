#ifndef Y3_CLUSTER_MAKE_INTERP_1D_CUH
#define Y3_CLUSTER_MAKE_INTERP_1D_CUH

#include "cosmosis/datablock/datablock.hh"
#include "common/cuda/Interp1D.cuh"

#include <vector>

using doubles = std::vector<double>;

inline quad::Interp1D
make_Interp1D(cosmosis::DataBlock& cfg,
              char const* section,
              char const* x_axis_name,
              char const* y_axis_name)
{
  std::vector<double> xs = cfg.view<doubles>(section, x_axis_name);
  std::vector<double> ys = cfg.view<doubles>(section, y_axis_name);
  assert(xs.size() == ys.size());
  return {xs.data(), ys.data(), xs.size()};
}

#endif
