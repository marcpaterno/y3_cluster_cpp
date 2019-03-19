#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "test/interp_1d.hh"
#include "test/primitives.hh"

namespace y3_cluster {

  class sigma_crit_inv {
  public:
    sigma_crit_inv(std::shared_ptr<Interp1D const> da, double c, double G)
      :_da(da), _c(c), _G(G)
    {}

    double
    operator()(double zt, zs) const
    {
      double _sig_crit_inv=0;
      if (zs > zt) {
        double const da_zt = _da->eval(zt), // da_z needs to be in Mpc
                     da_zs = _da->eval(zs), // da_z needs to be in Mpc
                     da_zt_zs = da_zs - (1.0+zt)/(1.0+zs) *da_zt; // da_z needs to be in Mpc
        _sig_crit_inv = 4.0*pi()*_G/_c/_c * da_zt * da_zt_zs / da_zs;
      }
      return _sig_crit_inv;
    }

  private:
    std::shared_ptr<Interp1D const> _da;
    double _c;
    double _G;
  };
}
