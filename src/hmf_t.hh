#ifndef Y3_CLUSTER_CPP_HMF_T_HH
#define Y3_CLUSTER_CPP_HMF_T_HH

#include <datablock_reader.hh>
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
          get_datablock<doubles>(sample, "HMF_params", "xs"),
          get_datablock<doubles>(sample, "HMF_params", "ys"),
          get_datablock<doubles>(sample, "HMF_params", "zs")))
      , _s(get_datablock<double>(sample, "MHF_params", "s"))
      , _q(get_datablock<double>(sample, "MHF_params", "q"))
    {}

    double
    operator()(double lnM, double zt) const
    {
      return _nmz->eval(lnM, zt) * (_s * (lnM - 37.5) + _q);
    }

  private:
    std::shared_ptr<Interp2D const> _nmz;
    double _s;
    double _q;
  };
}

#endif
