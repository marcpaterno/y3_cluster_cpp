#ifndef Y3_CLUSTER_CPP_HMF_T_HH
#define Y3_CLUSTER_CPP_HMF_T_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "interp_2d.hh"

#include <memory>

namespace y3_cluster {

  class HMF_t {
  public:
    HMF_t(std::shared_ptr<Interp2D const> nmz, double s, double q)
      : _nmz(nmz), _s(s), _q(q)
    {}

    using doubles = std::vector<double>;

    explicit HMF_t(cosmosis::DataBlock& sample)
      : _nmz(std::make_shared<Interp2D const>(
          sample.view<doubles>("HMF_params", "xs"),
          sample.view<doubles>("HMF_params", "ys"),
          sample.view<doubles>("HMF_params", "zs")))
      , _s(sample.view<double>("MHF_params", "s"))
      , _q(sample.view<double>("MHF_params", "q"))
    {}

    double
    operator()(double lnM, double zt) const
    {
      return _nmz->eval( lnM, zt) * (_s * ( lnM *0.4342944819 - 13.8) + _q);
      //0.4342944819 is log(e)
    }

  private:
    std::shared_ptr<Interp2D const> _nmz;
    double _s;
    double _q;
  };
}

#endif
