#ifndef Y3_CLUSTER_INT_LC_LT_DES_F_CUH
#define Y3_CLUSTER_INT_LC_LT_DES_F_CUH

#include "utils/cuda_interp_1d.cuh"
#include "cosmosis/datablock/datablock.hh"

// This is an error function fitted from the interpolation tables
// the parameter values were fitted follow a second order polynomial 
// this model is accurate to 1e-2

// to re-do the model fit go to ./y3_buzzard/fit_interpolator_int_lc_lt.py

namespace y3_cuda {
  // Define the interpolator values for the LC_LT_DES_f model
  // lbd_min, lbd_max and sigma
  // the function is the (erf[(xmax-x)/base]-erf[(xmin-x)/base])/2
  // where base = sqrt(2)*sigma
  // pivot redshift
  double zpivot = 0.35;
  std::array<std::array<double, 3>, 4> lminV = {{{19.12589544, -0.67211418, -3.28571528},
                                                 {27.33119808, -1.96743550, -6.08255365},
                                                 {39.09022658, -4.66476234, -9.55930173},
                                                 {49.06873402, -7.27724173, -10.92494645}}};

  std::array<std::array<double, 3>, 4> lmaxV = {{{26.50098912, -1.95380207, -5.06509283},
                                                {38.12144514, -4.15451144, -8.49090243},
                                                {50.00197827, -7.02875751, -12.40561456},
                                                {149.20116557, -26.90786363, -8.16398602}}};

  std::array<std::array<double, 3>, 4> sigmaV = {{{4.54634467, 1.45736381, -0.72576264},
                                                  {6.13291022, 1.47749438, -2.74077706},
                                                  {8.20122404, 2.42418783, -5.62164228},
                                                  {11.32534271, 4.66941640, -16.15420887}}};

  class INT_LC_LT_DES_f{
  public:
    INT_LC_LT_DES_f(): 
      _lminV(lminV), _lmaxV(lmaxV), _sigmaV(sigmaV), _zpivot(zpivot)
      { }

    explicit INT_LC_LT_DES_f(cosmosis::DataBlock const&): INT_LC_LT_DES_f() { }

    __device__ __host__ double
    operator()(int lo, double lt, double zt) const
    {
      // for a given lambda obs it finds the index of the lambda bin
      int index = (lo >= 60) * 3 + (lo >= 45 && lo < 60) * 2 + (lo >= 30 && lo < 45) * 1;

      double zp = zt-_zpivot;
      double l1 = _lminV[index][0] + _lminV[index][1]*zp + _lminV[index][2]*zp*zp;
      double l2 = _lmaxV[index][0] + _lmaxV[index][1]*zp + _lmaxV[index][2]*zp*zp;
      double s = _sigmaV[index][0] + _sigmaV[index][1]*zp + _sigmaV[index][2]*zp*zp;

      double base = sqrt(2.0)*s;
      double xu = l2-lt;
      double xl = l1-lt;

      double const val = ( erf(xu/base) - erf(xl/base) )/2.0;
      return val;
    }
  private:
    double _zpivot;
    std::array<std::array<double, 3>, 4> _lminV;
    std::array<std::array<double, 3>, 4> _lmaxV;
    std::array<std::array<double, 3>, 4> _sigmaV;
  };
}

#endif
