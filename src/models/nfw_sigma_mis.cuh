// Off-Centered Sigma NFW Profile
// Uses an interpolation table (look at data/nfw_off_center/)
// Assumes that datablock has rho_c, concetration
#ifndef Y3_CLUSTER_NFW_SIGMA_MIS
#define Y3_CLUSTER_NFW_SIGMA_MIS

#include <algorithm>
#include <iostream>
#include <math.h>
#include <string>

#include "cosmosis/datablock/datablock.hh"

#include "utils/make_interp_1d.hh"
#include "utils/cuda_interp_2d.cuh"


namespace y3_cuda {
  // Default concentration value
  double const CONC = 4.0;

  // Critical density in Msun/Mpc^3
  double const RHOC = 2.77533742639e+11;

  // selects the miscentering kernel ('single','gamma')
  std::string const GAMMA = "gamma";


  // Helper functions to construct filenames needed to read the interpolation
  // table information.
  inline std::string logx_file(std::string const& kernel)
  {
    return fmt::format("nfw_off_center/table_1000_1e-02_1e+04_{}_logx.txt", kernel);
  }

  inline std::string logxmis_file(std::string const& kernel)
  {
    return fmt::format("nfw_off_center/table_1000_1e-02_1e+04_{}_logxmis.txt", kernel);
  }

  inline std::string log_sigma_file(std::string const& kernel)
  {
    return fmt::format("nfw_off_center/table_1000_1e-02_1e+04_log_sigma_{}.txt", kernel);
  }


  class NFW_SIGMA_MIS {

    public:
    NFW_SIGMA_MIS(double c, double rhoc, std::string const& kernel)
      : _c(c),
        _rhoc(rhoc),
        _nfwProfile(read_vector(logx_file(kernel)),
                    read_vector(logxmis_file(kernel)),
                    read_vector(log_sigma_file(kernel)))
    { }

    NFW_SIGMA_MIS()
    : _c(CONC),
      _rhoc(RHOC),
      _nfwProfile(read_vector(logx_file(GAMMA)),
                  read_vector(logxmis_file(GAMMA)),
                  read_vector(log_sigma_file(GAMMA)))

    { }

    // In case, we envolve the NFW profile with redshift
    // TODO: Implement Mass-Concentration Relation
    // TODO: Implement different operator in case of rhocz(zt)
    // Ask Marc How to make _c and _rhoc be functional forms in any case
    NFW_SIGMA_MIS(cosmosis::DataBlock& sample)
    : _c(y3_cluster::make_Interp1D(sample,"correlationFunction","lnM","concentration").clamp(14.0))
    , _rhoc(y3_cluster::make_Interp1D(sample,"correlationFunction","z","rhoc").clamp(0.0))
    , _nfwProfile(read_vector(logx_file(GAMMA)),
                  read_vector(logxmis_file(GAMMA)),
                  read_vector(log_sigma_file(GAMMA)))
    { }

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

      double const log_unfw = _nfwProfile.clamp(std::log(x), std::log(xmis));
      
      double const nfw = norm * std::exp(log_unfw);
      return nfw;
    }

  private:
    double const _c;
    double const _rhoc;
    gpu_support::Interp2D _nfwProfile;
  };
}
#endif
