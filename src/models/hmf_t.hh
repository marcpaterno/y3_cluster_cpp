#ifndef Y3_CLUSTER_CPP_HMF_T_HH
#define Y3_CLUSTER_CPP_HMF_T_HH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "cosmosis/datablock/section_names.h"
#include "utils/interp_2d.hh"
#include "utils/read_vector.hh"

#include <memory>

namespace y3_cluster {

  namespace {
    inline std::vector<double>
    _adjust_to_log(cosmosis::DataBlock& db, const std::vector<double>& masses)
    {
      std::vector<double> output = masses;
      double const omega_m = db.view<double>("cosmological_parameters", "omega_M");
      double const omega_mu = db.view<double>("cosmological_parameters", "omega_nu");
      for (auto& x : output) x = std::log(x * (omega_m - omega_mu));
      return output;
    }

  }

  class HMF_t {
  public:
    HMF_t(Interp2D const& nmz, double s, double q)
      : _nmz(nmz), _s(s), _q(q)
    {}

    using doubles = std::vector<double>;

    explicit HMF_t(cosmosis::DataBlock& sample)
      : _nmz(_adjust_to_log(sample, sample.view<doubles>("mass_function", "m_h"))
            ,sample.view<doubles>("mass_function", "z")
            ,sample.view<cosmosis::ndarray<double>>("mass_function", "dndlnmh"))
      , _s(sample.view<double>("cluster_abundance", "hmf_s"))
      , _q(sample.view<double>("cluster_abundance", "hmf_q"))
    {}

    double
    operator()(double lnM, double zt) const
    {
      return _nmz.clamp(lnM, zt) *
             (_s * (lnM * 0.4342944819 - 13.8124426028) + _q);
      // 0.4342944819 is log(e)
    }

  private:
    Interp2D _nmz;
    double _s;
    double _q;
  };
}

#endif
