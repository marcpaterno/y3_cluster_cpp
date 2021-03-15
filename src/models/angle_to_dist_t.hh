#ifndef Y3_CLUSTER_ANGLE_TO_DIST_HH
#define Y3_CLUSTER_ANGLE_TO_DIST_HH

#include <vector>

#include "cosmosis/datablock/datablock.hh"
#include "utils/datablock_reader.hh"
#include "utils/interp_1d.hh"
#include "utils/primitives.hh"

namespace y3_cluster {

  class ANGLE_TO_DIST_t {
  public:
    ANGLE_TO_DIST_t(Interp1D da, double h0)
      : _da(std::move(da)), _h(h0)
    {}

    explicit ANGLE_TO_DIST_t(cosmosis::DataBlock& sample)
      : _da(get_datablock<std::vector<double>>(sample, "distances", "z"),
           get_datablock<std::vector<double>>(sample, "distances", "d_a"))
      , _h(get_datablock<double>(sample, "cosmological_parameters", "h0"))
    {}

    double
    operator()(double theta, double zt) const
    {
      double const da_zt = _da(zt) * _h; // da_zt in in Mpc/h
      double const dist =
        theta / 60.0 / 180.0 * pi() * da_zt * (1 + zt); // theta in arcminutes
      return dist;
    }

  private:
    Interp1D _da;
    double _h;
  };
}

#endif
