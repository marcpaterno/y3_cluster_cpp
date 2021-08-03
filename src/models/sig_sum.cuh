#ifndef Y3_CLUSTER_DEL_SIG_TOM_CUH
#define Y3_CLUSTER_DEL_SIG_TOM_CUH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "models/ez.cuh"
#include "cudaPagani/quad/GPUquad/Interp2D.cuh"
#include "utils/make_interp_2d.cuh"
#include "utils/primitives.cuh"

  class SIG_SUM {
  private:
    quad::Interp2D _sigma1;
    quad::Interp2D _sigma2;
    quad::Interp2D _bias;

  public:
    SIG_SUM(quad::Interp2D const& sigma1,
            quad::Interp2D const& sigma2,
            quad::Interp2D const& bias)
      : _sigma1(sigma1), _sigma2(sigma2), _bias(bias)
    {}

    using doubles = std::vector<double>;

    explicit SIG_SUM(cosmosis::DataBlock& sample)
      : _sigma1(make_Interp2D(sample,
                              "deltasigma",
                              "r_sigma_deltasigma",
                              "lnM",
                              "sigma_1"))
      , _sigma2(make_Interp2D(sample,
                              "deltasigma",
                              "r_sigma_deltasigma",
                              "matter_power_lin",
                              "z",
                              "deltasigma",
                              "sigma_2"))
      , _bias(make_Interp2D(sample,
                            "matter_power_lin",
                            "z",
                            "deltasigma",
                            "lnM",
                            "deltasigma",
                            "bias"))
    {}

    __device__
    double
    operator()(double r, double lnM, double zt) const
    /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
    {
      double _sig_1 = _sigma1.clamp(r, lnM);
      double _sig_2 = _bias.clamp(zt, lnM) * _sigma2.clamp(r, zt);
      // TODO: h factor?
      // if (del_sig_1 >= del_sig_2) {
      // return (1.+zt)*(1.+zt)*(1.+zt)*(_sig_1+_sig_2);
      return (_sig_1 + _sig_2);
      /*} else {
        return 1e12*(1.+zt)*(1.+zt)*(1.+zt)*del_sig_2;
      } */
    }
  };

#endif
