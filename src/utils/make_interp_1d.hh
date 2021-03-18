#ifndef Y3_CLUSTER_MAKE_INTERP_1D_HH
#define Y3_CLUSTER_MAKE_INTERP_1D_HH

#include "utils/interp_1d.hh"
#include "cosmosis/datablock/datablock.hh"

#include <vector>

namespace y3_cluster {

  using doubles = std::vector<double>;
  
  inline
  Interp1D make_Interp1D(cosmosis::DataBlock& cfg,
      char const* section,
      char const* x_axis_name,
      char const* y_axis_name) {
    return Interp1D(cfg.view<doubles>(section, x_axis_name),
                    cfg.view<doubles>(section, y_axis_name));
  }
}


#endif
