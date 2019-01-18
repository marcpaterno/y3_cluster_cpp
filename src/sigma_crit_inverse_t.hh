#include <exception>
#include <memory>
#include <vector>

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "datablock_reader.hh"
#include "interp_1d.hh"
#include "primitives.hh"

namespace y3_cluster {

  class sigma_crit_inv {
  public:
    sigma_crit_inv(std::shared_ptr<Interp1D const> da)
      : _da(da)
    {}

    explicit sigma_crit_inv(cosmosis::DataBlock& sample)
      : sigma_crit_inv(std::make_shared<y3_cluster::Interp1D const>(
                     get_datablock<std::vector<double>>(sample, "distances", "z"),
                     get_datablock<std::vector<double>>(sample, "distances", "d_a")))
    {}

    sigma_crit_inv(const sigma_crit_inv&) = default;
    sigma_crit_inv(sigma_crit_inv&&) = default;

    double
    operator()(double zt, double zs) const
    {
      double _sig_crit_inv = 0;
      if (zs > zt) {
      // Units in Mpc, M_sol, s
        double const da_zt = _da->eval(zt), // da_z needs to be in Mpc
                     da_zs = _da->eval(zs), // da_z needs to be in Mpc
		     // TODO: Check this - is likely wrong!
                     da_zt_zs = da_zs - (1.0+zt)/(1.0+zs) *da_zt; // da_z needs to be in Mpc
        _sig_crit_inv = 4.0*pi()*G()/c()/c() * da_zt * da_zt_zs / da_zs;
      }
      return _sig_crit_inv;
    }

  private:
    std::shared_ptr<Interp1D const> _da;
    double _c;
    double _G;
  };
}
