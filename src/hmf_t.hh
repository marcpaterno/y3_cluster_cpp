#ifndef Y3_CLUSTER_CPP_HMF_T_HH
#define Y3_CLUSTER_CPP_HMF_T_HH

#include <datablock_reader.hh>
#include "interp_2d.hh"
#include <read_vector.hh>
#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "/cosmosis/cosmosis/datablock/section_names.h"

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
                  get_datablock<doubles>(sample, "mass_function", "m_h"),
                  get_datablock<doubles>(sample, "mass_function", "z"),
                  read_vector("dndlnmh.txt", [](double x) { return x; })))
                  //get_datablock<std::vector<doubles>>(sample, "mass_function", "dndlnmh")))
      , _s(get_datablock<double>(sample, "gamma_t", "hmf_s"))
      , _q(get_datablock<double>(sample, "gamma_t", "hmf_q"))
    {}

    double
    operator()(double lnM, double zt) const
    {
      return _nmz->eval( lnM, zt) * (_s * ( lnM *0.4342944819 - 13.8124426028) + _q);
      //0.4342944819 is log(e)
    }

  private:
    std::shared_ptr<Interp2D const> _nmz;
    double _s;
    double _q;
  };
}

#endif
