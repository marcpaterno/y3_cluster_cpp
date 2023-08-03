// Off-Centered Sigma NFW Profile
// Uses an interpolation table (look at data/nfw_off_center/)
// Assumes that datablock has rho_c, concetration
#ifndef Y3_CLUSTER_NFW_DSIGMA_MIS
#define Y3_CLUSTER_NFW_DSIGMA_MIS

#include <algorithm>
#include <iostream>
#include <math.h>
#include <string>
#include "cosmosis/datablock/ndarray.hh"
#include "common/cuda/Interp2D.cuh"
#include "utils/make_interp_1d.cuh"
#include "utils/primitives.hh"
#include "utils/read_vector.hh"

namespace y3_cuda {
  // Default concentration value
  double const CONC = 4.0;

  // Critical density in Msun/Mpc^3
  double const RHOC = 2.77533742639e+11;

  // selects the miscentering kernel ('single','gamma')
  std::string const GAMMA = "gamma";
  
  class NFW_DSIGMA_MIS {

    public:
    NFW_DSIGMA_MIS(double c, double rhoc, std::string const& kernel)
      : _c(c),
        _rhoc(rhoc)
    {
      std::string xfile = "nfw_off_center/offset_nfw_table_500_1e-02_1e+04_" + kernel + "_logx.txt";
      std::string yfile = "nfw_off_center/offset_nfw_table_500_1e-02_1e+04_" + kernel + "_logxmis.txt";
      std::string zfile = "nfw_off_center/offset_nfw_table_500_1e-02_1e+04_log_deltasigma_" + kernel + ".txt";
		
      auto const xs = read_vector(xfile);
      auto const ys = read_vector(yfile);
      auto const zs = read_vector(zfile);
      quad::Interp2D temp(xs, ys, zs);

      _nfwProfile = temp;
    }

    explicit NFW_DSIGMA_MIS()
    : _c(CONC),
      _rhoc(RHOC)
    {
      std::string xfile = "nfw_off_center/offset_nfw_table_500_1e-02_1e+04_"+ GAMMA + "_logx.txt";
      std::string yfile = "nfw_off_center/offset_nfw_table_500_1e-02_1e+04_"+ GAMMA + "_logxmis.txt";
      std::string zfile = "nfw_off_center/offset_nfw_table_500_1e-02_1e+04_log_deltasigma_" + GAMMA + ".txt";
		
      auto const xs = read_vector(xfile);
      auto const ys = read_vector(yfile);
      auto const zs = read_vector(zfile);

      quad::Interp2D temp(xs, ys, zs);
      _nfwProfile = temp;
    }

    // In case, we envolve the NFW profile with redshift
    // TODO: Implement Mass-Concentration Relation
    // TODO: Implement different operator in case of rhocz(zt)
    // Ask Marc How to make _c and _rhoc be functional forms in any case
    // NFW_DSIGMA_MIS(cosmosis::DataBlock& sample)
    // : _c(make_Interp1D(sample,"correlationFunction","lnM","concentration").clamp(14.0))
    // , _rhoc(make_Interp1D(sample,"correlationFunction","z","rhoc").clamp(0.0))
    // {
    //   std::string xfile = "nfw_off_center/offset_nfw_table_500_1e-02_1e+04_"+ GAMMA + "_logx.txt";
    //   std::string yfile = "nfw_off_center/offset_nfw_table_500_1e-02_1e+04_"+ GAMMA + "_logxmis.txt";
    //   std::string zfile = "nfw_off_center/offset_nfw_table_500_1e-02_1e+04_log_deltasigma_" + GAMMA + ".txt";
		
    //   auto const xs = read_vector(xfile);
    //   auto const ys = read_vector(yfile);
    //   auto const zs = read_vector(zfile);
    //   quad::Interp2D temp = quad::Interp2D(xs, ys, zs);
    //   _nfwProfile = temp;
    // }

    __device__ __host__ double
    operator()(double r, double rmis, double lnM) const 
    {
      double const delta_c = 200.0 * _c * _c * _c / 3.0 * (std::log(1.0 + _c) - _c / (1.0 + _c));
      double const rho_crit = _rhoc;
      double const r_200 = std::cbrt(3.0 * std::exp(lnM) / (800.0 * M_PI * rho_crit));
      double const r_s = r_200 / _c;
      double const norm = r_s * delta_c * rho_crit;

      double const x = r / r_s;
      double const xmis = rmis / r_s;

      double log_unfw = 0;
      
#ifdef __CUDA_ARCH__
        log_unfw = _nfwProfile.clamp(std::log(x), std::log(xmis));
#else
        // THIS IS DUMMY CODE. IT WILL BE REPLACED.
        log_unfw = 0.;
#endif
      
      double const nfw = norm * std::exp(log_unfw);
      return nfw;
    }

    size_t
    get_device_mem_footprint()
    {
      return _nfwProfile.get_device_mem_footprint();
    }


  private:
    double const _c;
    double const _rhoc;
    quad::Interp2D _nfwProfile;
  };
}
#endif
