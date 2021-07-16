#ifndef DV_DO_DZ_T_CUH
#define DV_DO_DZ_T_CUH

#include "cudaPagani/quad/Interp1D.cuh"
#include "models/ez.cuh"

class dv_do_dz_t {
public:
  dv_do_dz_t() = default;

  dv_do_dz_t(quad::Interp1D* da, ez ezt, double h)
    : _da(da), _ezt(ezt), _h(h)
  {}

  using doubles = std::vector<double>;

  __device__ double
  operator()(double zt) const
  {
    double const da_z = _da->eval(zt); // da_z needs to be in Mpc
    // Units: (Mpc/h)^3
    // 2997.92 is Hubble distance, c/H_0
    return 2997.92 * (1.0 + zt) * (1.0 + zt) * da_z * _h * da_z * _h / _ezt(zt);
  }

private:
  quad::Interp1D* _da;
  ez _ezt;
  double _h;
};

#endif
