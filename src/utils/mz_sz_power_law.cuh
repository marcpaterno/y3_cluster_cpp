#ifndef Y3_CLUSTER_MZ_POWER_LAW_CUH
#define Y3_CLUSTER_MZ_POWER_LAW_CUH

// mz_sz_power_law represents a commonly-used power law relationship used
// in SPT related works (Benson et al. (2013))
//       log(A) + B log(M/Mp) + C log(E(z)/E(z0))
// with A, B and C being constants set in the construction of the power law
// object.
// And Mp=3 10^14 and z0=0.6

class mz_sz_power_law {
public:
  mz_sz_power_law(double A, double B, double C, double om) noexcept;

  // The function call operator evaluates the power law at the given
  // values of lnM and z. Note that the first parameter is not mass,
  // but ln(mass).
  __device__ __host__ double operator()(double lnM, double z) const noexcept;

private:
  double const log_A_;
  double const B_;
  double const C_;
  double const omega_m_;
};

inline mz_sz_power_law::mz_sz_power_law(double A, double B, double C, double om) noexcept
  : log_A_(log(A)), B_(B), C_(C), omega_m_(om)
{}

inline __device__ __host__ double
mz_sz_power_law::operator()(double lnM, double z) const noexcept
{
  // pivot values
  double const Mp = 3.0e14
  double const z0 = 0.6

  // Redshift evolution term
  double const omega_l_ = 1-omega_m_
  double const Ez0 = sqrt(omega_m_ * (z0+1.) * (z0+1.) * (z0+1.) + omega_l_)
  double const Ez = sqrt(omega_m_ * (z+1.) * (z+1.) * (z+1.) + omega_l_)

  // Mass-Observarble Relation
  double const log_res =  log_A_ + B_ * (lnM-log(Mp)) + C_ * log(Ez/Ez0);
  return exp(log_res);
}

#endif
