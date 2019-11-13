#ifndef Y3_CLUSTER_CPP_MOR_1D_HH
#define Y3_CLUSTER_CPP_MOR_1D_HH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "utils/datablock_reader.hh"
#include "utils/interp_2d.hh"
#include "utils/primitives.hh"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

namespace y3_cluster {

  using doubles = std::vector<double>;

  class MOR_1D {
  public:
    MOR_1D(double Mmin,
           double M1,
           double alpha,
           double epsilon,
           double sigintr,
           double zplam)
      : Mmin_(Mmin)
      , M1_(M1)
      , alpha_(alpha)
      , epsilon_(epsilon)
      , sigintr_(sigintr)
      , zplam_(zplam)
    {}

    MOR_1D(cosmosis::DataBlock& sample)
      : Mmin_(get_datablock<double>(sample, "abundance_integral", "Mmin"))
      , M1_(get_datablock<double>(sample, "abundance_integral", "M1"))
      , alpha_(get_datablock<double>(sample, "abundance_integral", "alpha"))
      , epsilon_(get_datablock<double>(sample, "abundance_integral", "epsilon"))
      , sigintr_(get_datablock<double>(sample, "abundance_integral", "sig0lam"))
      , zplam_(get_datablock<double>(sample, "abundance_integral", "zplam"))
    {}

    double
    operator()(double lamtrue, double ztrue, double lnM200m) const
    {
      double const mean_lamsat =
        pow((std::exp(lnM200m) - Mmin_) / (M1_ - Mmin_), alpha_) *
        pow((1.0 + ztrue) / (1.0 + zplam_), epsilon_);
      double const mean_lamtrue = 1.0 + mean_lamsat;
      double const var_lam =
        mean_lamsat + mean_lamsat * mean_lamsat * sigintr_ * sigintr_;

      return y3_cluster::gaussian(lamtrue, mean_lamtrue, std::sqrt(var_lam));
    }

  private:
    double const Mmin_, M1_, alpha_, epsilon_, sigintr_, zplam_;
  };
}

#endif
