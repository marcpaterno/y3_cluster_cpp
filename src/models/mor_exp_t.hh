#ifndef Y3_CLUSTER_CPP_MOR_T_HH
#define Y3_CLUSTER_CPP_MOR_T_HH

#include "utils/datablock_reader.hh"
#include "utils/mz_power_law.hh"
#include "utils/primitives.hh"

#include <cmath>

namespace y3_cluster {

  class MOR_EXP_t {
  public:
    MOR_EXP_t(double A,
              double B,
              double C,
              double sigma_i,
              double epsilon,
              double skew,
              double z_pivot)
      : _A(A)
      , _B(B)
      , _C(C)
      , _sigma_intr(sigma_i)
      , _epsilon(epsilon)
      , _z_pivot(z_pivot)
      , _skew(skew)
    {}

    explicit MOR_EXP_t(cosmosis::DataBlock& sample)
      : _A(get_datablock<double>(sample, "cluster_abundance", "mor_A"))
      , _B(get_datablock<double>(sample, "cluster_abundance", "mor_B"))
      , _C(get_datablock<double>(sample, "cluster_abundance", "mor_alpha"))
      , _sigma_intr(
          get_datablock<double>(sample, "cluster_abundance", "mor_sigma"))
      , _epsilon(
          get_datablock<double>(sample, "cluster_abundance", "mor_epsilon"))
      , _z_pivot(
          get_datablock<double>(sample, "cluster_abundance", "z_mor_pivot"))
      , _skew(get_datablock<double>(sample, "cluster_abundance", "mor_skew"))
    {}

    double
    operator()(double lt, double lnM, double zt) const
    {
      double ltm = pow((std::exp(lnM) - _A) / (_B - _A), _C) *
                   pow((1.0 + zt) / (1.0 + _z_pivot), _epsilon);
      double const x = lt - ltm;
      double const erfarg = -1.0 * _skw * (x) / (std::sqrt(2.) * _sigma);
      double const erfterm = std::erfc(erfarg);
      return y3_cluster::gaussian(x, 0.0, _sigma_intr) * erfterm;
    }

  private:
    double _A;
    double _B;
    double _C;
    double _sigma_intr;
    double _skw;
    double _epsilon;
    double _z_pivot
  };
}

#endif
