#ifndef Y3_CLUSTER_LC_LT_T_HH
#define Y3_CLUSTER_LC_LT_T_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include <cmath>

namespace y3_cluster {
  class LC_LT_t {
  public:
    LC_LT_t(double tau, double mu, double sigma, double fmsk, double fprj)
      : _tau(tau), _mu(mu), _sigma(sigma), _fmsk(fmsk), _fprj(fprj)
    {}

    explicit LC_LT_t(cosmosis::DataBlock& sample)
    {
      sample.get_val<double>("LC_LT_params", "tau", _tau);
      sample.get_val<double>("LC_LT_params", "mu", _mu);
      sample.get_val<double>("LC_LT_params", "sigma", _sigma);
      sample.get_val<double>("LC_LT_params", "fmsk", _fmsk);
      sample.get_val<double>("LC_LT_params", "fprj", _fprj);
    }

    double
    operator()(double lc, double lt, double /* zt */) const
    {
      return (1.0 - _fmsk) * (1.0 - _fprj) * y3_cluster::invsqrt2pi() *
               std::exp(-(lc - _mu) * (lc - _mu) / (2.0 * _sigma * _sigma)) /
               _sigma +
             0.5 * ((1.0 - _fmsk) * _fprj * _tau + _fmsk * _fprj / lt) *
               std::exp(_tau * (2.0 * _mu + _tau * _sigma * _sigma - 2.0 * lc) /
                        2.0) *
               std::erfc((_mu + _tau * _sigma * _sigma - lc) /
                         (std::sqrt(2.0) * _sigma)) +
             _fmsk *
               (std::erfc((_mu - lc - lt) / (std::sqrt(2.0) * _sigma)) -
                std::erfc((_mu - lc) / (std::sqrt(2.0) * _sigma))) /
               (2.0 * lt) -
             _fmsk * _fprj *
               (std::exp(-_tau * lt) *
                std::exp(_tau *
                         (2.0 * _mu + _tau * _sigma * _sigma - 2.0 * lc) /
                         2.0) *
                std::erfc((_mu + _tau * _sigma * _sigma - lc - lt) /
                          (std::sqrt(2.0) * _sigma))) /
               (2.0 * lt);
    }

  private:
    double _tau;
    double _mu;
    double _sigma;
    double _fmsk;
    double _fprj;
  };
}

#endif