#ifndef Y3_CLUSTER_LN_MEZ_POWER_LAW_HH
#define Y3_CLUSTER_LN_MEZ_POWER_LAW_HH

#include "ez.hh"
#include "utils/datablock_reader.hh"

#include <cmath>

namespace y3_cluster {
  class ln_mez_power_law {
  public:
    ln_mez_power_law(double A,
                     double B,
                     double C,
                     double lnMp,
                     double zp,
                     cosmosis::DataBlock& sample)
      : log_A_(std::log(A))
      , B_(B)
      , C_(C)
      , lnMp_(lnMp)
      , ez_(y3_cluster::EZ(sample))
      , ln_e_zp_(std::log(ez_(zp)))
    {}

    ln_mez_power_law(double A,
                     double B,
                     double C,
                     double lnMp,
                     double zp,
                     double omega_m,
                     double omega_l,
                     double omega_k)
      : log_A_(std::log(A))
      , B_(B)
      , C_(C)
      , lnMp_(lnMp)
      , ez_(omega_m, omega_l, omega_k)
      , ln_e_zp_(std::log(ez_(zp)))
    {}

    // ln( A * (m/mp)**B * (E(z)/E(zp))**C )
    double
    operator()(double lnM, double z) const
    {
      return log_A_ + B_ * (lnM - lnMp_) + C_ * (std::log(ez_(z)) - ln_e_zp_);
    }

  private:
    double const log_A_;
    double const B_;
    double const C_;
    double const lnMp_;
    y3_cluster::EZ ez_;
    double const ln_e_zp_;
  };
}
#endif
