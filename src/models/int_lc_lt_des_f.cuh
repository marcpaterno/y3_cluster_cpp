#ifndef Y3_CLUSTER_INT_LC_LT_DES_F_CUH
#define Y3_CLUSTER_INT_LC_LT_DES_F_CUH

#include "utils/cuda_interp_1d.cuh"
#include "cosmosis/datablock/datablock.hh"

// This is an error function fitted from the interpolation tables
// the parameter values were fitted follow a second order polynomial 
// this model is accurate to 5e-3

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

  std::vector<double> lbdBins = {20, 30, 45, 60};
  std::vector<double> indBins = {0, 1, 2, 3};
  gpu_support::Interp1D interpInd(lbdBins, indBins);

  class INT_LC_LT_DES_f {
  public:
    // INT_LC_LT_DES_f(lbdBins, indBins): { }
    explicit INT_LC_LT_DES_f(cosmosis::DataBlock const&) { }

    __device__ __host__ double
    operator()(int lo, double lt, double zt) const
    {
      // Perform interpolation and convert to int
      double index_double = interpInd.clamp(lo);
      int index = static_cast<int>(index_double);

      // double val = 0;
      std::array<double, 3> lmin = lminV[index];
      std::array<double, 3> lmax = lmaxV[index];
      std::array<double, 3> sigma = sigmaV[index];
      
      double z = zt-zpivot;
      double l1 = lmin[0] + lmin[1]*z + lmin[2]*z*z;
      double l2 = lmax[0] + lmax[1]*z + lmax[2]*z*z;
      double s = sigma[0] + sigma[1]*z + sigma[2]*z*z;    

      double base = sqrt(2.0)*s;
      double xu = l2-lt;
      double xl = l1-lt;

      double const val = ( erf(xu/base) - erf(xl/base) )/2.0;
      return val;
    }
  };
}

#endif
