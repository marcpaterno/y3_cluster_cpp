#ifndef EZ_SQ_CUH
#define EZ_SQ_CUH

class ez_sq {
public:
  ez_sq() = default;

  ez_sq(double omega_m, double omega_l, double omega_k)
    : _omega_m(omega_m), _omega_l(omega_l), _omega_k(omega_k)
  {}

  __device__ double
  operator()(double z) const
  {
    // NOTE: this is valid only for \Lambda CDM cosmology, not wCDM
    double const zplus1 = 1.0 + z;
    return (_omega_m * zplus1 * zplus1 * zplus1 +
            _omega_k * zplus1 * zplus1 +
            _omega_l);
  }

private:
  double _omega_m = 0.0;
  double _omega_l = 0.0;
  double _omega_k = 0.0;
};

endif
