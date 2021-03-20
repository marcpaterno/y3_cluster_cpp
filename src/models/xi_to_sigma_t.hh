#ifndef Y3_CLUSTER_SIG_T_HH
#define Y3_CLUSTER_SIG_T_HH

#include "cosmosis/datablock/datablock.hh"
#include "utils/interp_2d.hh"
#include "utils/make_interp_2d.hh"

#include <cmath>
#include <cubacpp/cubacpp.hh>
#include <iostream>

namespace y3_cluster {
  class SIG_y3 {
  private:
    Interp2D _xi;
    double _om;

  public:
    SIG_y3(Interp2D const& xi, double om) : _xi(xi), _om(om) {}

    using doubles = std::vector<double>;

    explicit SIG_y3(doubles dim0,
                    doubles dim_R,
                    doubles xi_v,
                    cosmosis::DataBlock& sample)
      : _xi(dim_R, dim0, xi_v)
      , _om(sample.view<double>("cosmological_parameters", "omega_m"))
    {}

    double
    operator()(double R, double zt) const
    /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
    {

      // double rho_m = 0.277526157*_om;
      double rho_m = 0.277536627 * _om;
      cubacpp::QAG integrator(0, 150, GSL_INTEG_GAUSS61, 50000000);
      const auto res = integrator.integrate(
        [=](double chi) {
          double const r = std::hypot(R, chi);
          if (r > 150.) return 0.0;
          return  _xi(r, zt);
        },
        5e-3,
        1e-9);

      return res.value * rho_m * 2.0;
    }
  };
}

#endif
