#ifndef Y3_CLUSTER_CPP_HMF_T_HH
#define Y3_CLUSTER_CPP_HMF_T_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "test/interp_1d.hh"

namespace y3_cluster {

  class HMF_t {
  public:
    HMF_t(Interp1D const* nmz, double s, double q) : _nmz(nmz), _s(s), _q(q) {}

    explicit HMF_t(cosmosis::DataBlock& sample)
    // TODO: test once Interp2D implementation is finished
    //: _nmz([](cosmosis::DataBlock& x) {
    //  std::vector<double> const& xs, ys;
    //  cosmosis::ndarray<double> const& zs;
    //  x.get_val<std::vector<double>>("HMF_params", "xs", xs);
    //  x.get_val<std::vector<double>>("HMF_params", "ys", ys);
    //  x.get_val<cosmosis::ndarray<double>>("HMF_params", "zs", zs);
    //  return Interp2D{xs, ys, zs};
    //}(sample))
    {
      sample.get_val<double>("HMF_params", "s", _s);
      sample.get_val<double>("HMF_params", "q", _q);
    }

    double
    operator()(double lnM, double /*zt*/) const
    {
      // TODO: This is clearly worng!
      return _nmz->eval(lnM) * (_s * (lnM - 37.5) + _q);
    }

  private:
    // TODO: change to Interp2D once implementation is finished
    Interp1D const* _nmz;
    double _s;
    double _q;
  };
}

#endif