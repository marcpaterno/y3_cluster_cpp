// Off-Centered Sigma NFW Profile
// Uses an interpolation table (look at data/nfw_off_center/)
// Assumes that datablock has rho_c, concetration
#ifndef Y3_CLUSTER_NFW_SIGMA_MIS
#define Y3_CLUSTER_NFW_SIGMA_MIS

#include <iostream>
#include <math.h>
#include <string>
#include "cosmosis/datablock/ndarray.hh"
#include "common/cuda/Interp2D.cuh"
#include "utils/primitives.hh"

#include <algorithm>

namespace y3_cuda {
  class NFW_SIGMA_MIS {
    std::vector<y3_cluster::Interp2D> nfwProfile;

    // Default concentration value
    double const CONC = 4.0;

    // Critical density in Msun/Mpc^3
    double const RHOC = 2.77533742639e+11;

    // selects the miscentering kernel ('single','gamma')
    std::string GAMMA = "gamma";

  public:
    NFW_SIGMA_MIS(double c, double rhoc, std::string kernel)
      : _c(c), _rhoc(rhoc), _kernel(kernel)
    {
      std::string xfile = "data/nfw_off_center/table_1000_1e-02_1e+04_"+ std::string(_kernel) + "_logx.txt";
      std::string yfile = "data/nfw_off_center/table_1000_1e-02_1e+04_"+ std::string(_kernel) + "_logxmis.txt";
      std::string zfile = "data/nfw_off_center/table_1000_1e-02_1e+04_log_sigma_" + std::string(_kernel) + ".txt";
		
      auto const xs = read_vector(xfile);
      auto const ys = read_vector(yfile);
      auto const zs = read_vector(zfile);
      nfwProfile = y3_cluster::Interp2D(xs, ys, zs);
    }
    using doubles = std::vector<double>;

    explicit NFW_SIGMA_MIS()
    : _c(CONC), _rhoc(RHOC), _kernel(GAMMA)
    {
      std::string xfile = "data/nfw_off_center/table_1000_1e-02_1e+04_"+ std::string(_kernel) + "_logx.txt";
      std::string yfile = "data/nfw_off_center/table_1000_1e-02_1e+04_"+ std::string(_kernel) + "_logxmis.txt";
      std::string zfile = "data/nfw_off_center/table_1000_1e-02_1e+04_log_sigma_" + std::string(_kernel) + ".txt";
		
      auto const xs = read_vector(xfile);
      auto const ys = read_vector(yfile);
      auto const zs = read_vector(zfile);
      nfwProfile = y3_cluster::Interp2D(xs, ys, zs);
    }

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

      double const log_unfw = nfwProfile.clamp(std::log(x), std::log(xmis));
      double const nfw = norm * std::exp(log_unfw);
      return nfw;
    }

  private:
    std::string _kernel;
    double const _c;
    double const _rhoc;
  };
}
#endif

// /global/common/software/des/jesteves/y3_cluster_cpp/data/nfw_off_center
// 'offset_nfw_table_1000_1e-02_1e+04_log_sigma_single'

// Based on ChatGPT
void testNFW_SIGMA_MIS() {
  // Create an instance of NFW_SIGMA_MIS with test parameters
  double c = 4.0;
  double rhoc = 2.77533742639e+11;
  std::string kernel = "gamma";
  y3_cuda::NFW_SIGMA_MIS nfw(c, rhoc, kernel);

  // Define test inputs and expected output
  double r = 1.0;
  double rmis = 0.5;
  double lnM = 10.0;
  double expected = /* Expected output value */;

  // Evaluate the NFW_SIGMA_MIS operator and compare with expected output
  double result = nfw(r, rmis, lnM);
  assert(result == expected);

  // Print test result
  std::cout << "Test Passed!" << std::endl;
}

int main() {
  testNFW_SIGMA_MIS();
  return 0;
}
