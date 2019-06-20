#ifndef Y3_CLUSTER_ANGLE_TO_DIST_HH
#define Y3_CLUSTER_ANGLE_TO_DIST_HH

#include <exception>
#include <memory>
#include <vector>

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "utils/datablock_reader.hh"
#include "utils/interp_1d.hh"
#include "utils/primitives.hh"

namespace y3_cluster {

  class ANGLE_TO_DIST_t {
  public:
    ANGLE_TO_DIST_t(std::shared_ptr<Interp1D const> da, double h0)
      : _da(da)
      , _h(h0)
    {}

    explicit ANGLE_TO_DIST_t(cosmosis::DataBlock& sample)
      : _da(std::make_shared<y3_cluster::Interp1D const>(
                     get_datablock<std::vector<double>>(sample, "distances", "z"),
                     get_datablock<std::vector<double>>(sample, "distances", "d_a")))
	,_h(get_datablock<double>(sample, "cosmological_parameters", "h0"))
    {}

    double
    operator()(double theta, double zt) const
    {
      double _dist = 0;
      double const da_zt = _da->eval(zt)*_h; // da_zt in in Mpc/h
      _dist = theta/60.0/180.0*pi() * da_zt *(1+zt); //theta in arcminutes
      return _dist;
    }

  private:
    std::shared_ptr<Interp1D const> _da;
    double _h;
  };
}

#endif
