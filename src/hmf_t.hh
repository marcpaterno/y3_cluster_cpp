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
      : _nmz([sample] () mutable {
	// FIXME: Make this grab the hmf array from cosmosis!
	  auto identity = [](double x) { return x; };
	  std::cout << "Has gamma_t section: " << sample.has_section("gamma_t") << std::endl;
	  std::cout << "Has cosmo section: " << sample.has_section(COSMOLOGICAL_PARAMETERS_SECTION) << std::endl;
	  std::cout << "Has omega_m: " << sample.has_val("cosmological_parameters", "omega_m") << std::endl;
	  std::cout << sample.view<double>("cosmological_parameters", "omega_m") << std::endl;
	  double om = get_datablock<double>(sample, "cosmological_parameters", "omega_m");
	  //double om = 0.3;
	  auto log_omega_m = [om](double x) { return std::log(x*om); };
	  auto mh = read_vector("m_h.txt", log_omega_m);
	  auto const zz = read_vector("z.txt", identity);
	  auto const dndlnmh = read_vector("dndlnmh.txt", identity);
          return std::make_shared<Interp2D const>(mh, zz, dndlnmh);}())
          //get_datablock<doubles>(sample, "gamma_t", "hmf_xs"),
          //get_datablock<doubles>(sample, "gamma_t", "hmf_ys"),
          //get_datablock<doubles>(sample, "gamma_t", "hmf_zs")))
      , _s(get_datablock<double>(sample, "gamma_t", "hmf_s"))
      , _q(get_datablock<double>(sample, "gamma_t", "hmf_q"))
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
