#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "utils/interp_2d.hh"
#include "utils/make_interp_2d.hh"
#include "utils/primitives.hh"

#include <cubacpp/cubacpp.hh>

#include <cmath>

namespace y3_cluster {

  class xi_nl {
  public:
    xi_nl(Interp2D const& pkz) : _pkz(pkz) {}

    explicit xi_nl(cosmosis::DataBlock& sample)
      : _pkz(make_Interp2D(sample, "matter_power_nl", "k_h", "z", "p_k"))
    {}

    double
    operator()(double r, double z) const
    {

      cubacpp::QAG integrator(0.0001, 100.0, GSL_INTEG_GAUSS61, 5000000);
      const auto res = integrator.integrate(
        [=](double k) {
          double const x = k * r;
          double const w = sin(x) / x;
          return w * k * k * _pkz(k, z) / (2.0 * pi() * pi());
        },
        1e-3,
        1e-5);
      return res.value;
    }

  private:
    Interp2D _pkz;
  };
}