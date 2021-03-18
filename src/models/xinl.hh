#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "cosmosis/datablock/section_names.h"
#include "utils/datablock_reader.hh"
#include "utils/interp_2d.hh"
#include <math.h>
#include <memory>
using namespace y3_cluster;
#define PI 3.14159265
#include <cubacpp/cubacpp.hh>

class xi_nl {
public:
  xi_nl(std::shared_ptr<Interp2D const> pkz) : _pkz(pkz) {}

  explicit xi_nl(cosmosis::DataBlock& sample)
    : _pkz(std::make_shared<Interp2D const>(
        get_datablock<std::vector<double>>(sample, "matter_power_nl", "k_h"),
        get_datablock<std::vector<double>>(sample, "matter_power_nl", "z"),
        get_datablock<cosmosis::ndarray<double>>(sample,
                                                 "matter_power_nl",
                                                 "p_k")))
  {}

  double
  operator()(double r, double z) const
  {

    cubacpp::QAG integrator(0.0001, 100.0, GSL_INTEG_GAUSS61, 5000000);
    const auto res = integrator.integrate(
      [=](double k) {
        double x = k * r;
        double w = sin(x) / x;
        double rest = w * k * k * _pkz->eval(k, z) / (2.0 * PI * PI);
        return rest;
      },
      1e-3,
      1e-5);
    return res.value;
  }

private:
  std::shared_ptr<Interp2D const> _pkz;
};
