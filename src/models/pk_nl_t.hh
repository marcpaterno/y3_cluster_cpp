#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "cosmosis/datablock/section_names.h"

#include "cosmosis/datablock/datablock.hh"
#include "utils/interp_2d.hh"
#include "utils/make_interp_2d.hh"

class pk_nl {
public:
  explicit pk_nl(y3_cluster::Interp2D const& pkz) : _pkz(pkz) {}

  explicit pk_nl(cosmosis::DataBlock& sample)
    : _pkz(y3_cluster::make_Interp2D(sample, "matter_power_nl", "k_h", "z", "p_k"))
  {}

  double
  operator()(double k, double z) const
  {

    return _pkz.clamp(k, z);
  }

private:
  y3_cluster::Interp2D _pkz;
};
