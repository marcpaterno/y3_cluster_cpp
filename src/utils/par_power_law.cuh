#ifndef Y3_CLUSTER_PAR_POWER_LAW_CUH
#define Y3_CLUSTER_PAR_POWER_LAW_CUH

#include <cmath> // For log and exp functions

namespace y3_cuda {

class par_power_law {
  double const ZPIV = 0.35;
  double const XPIV = 1.0;

public:
  // Constructors
  par_power_law(double A, double B, double C) noexcept
    : log_A_(std::log(A)), B_(B), C_(C), zp_(ZPIV), xp_(XPIV) {}

  par_power_law(double A, double B, double C, double zp, double xp) noexcept
    : log_A_(std::log(A)), B_(B), C_(C), zp_(zp), xp_(xp) {}

  // Function call operator
  __device__ __host__ double operator()(double lnM, double z) const noexcept {
    double const log_res = B_ * std::log(std::exp(lnM)/xp_) + C_ * std::log((1+z)/(1+zp_)) + log_A_;
    return std::exp(log_res);
  }

private:
  double const log_A_;
  double const B_;
  double const C_;
  double const zp_;
  double const xp_;
};

} // namespace y3_cuda

#endif
