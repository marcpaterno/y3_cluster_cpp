#ifndef ROFFSET_T_CUH
#define ROFFSET_T_CUH

class roffset_t {
public:
  roffset_t() = default;

  explicit roffset_t(double tau) : _tau(tau) {}

  __device__ double
  operator()(double x) const
  {
    // eq. 36
    return x / _tau / _tau * exp(-x / _tau);
  }

private:
  double _tau = 0.0;
};


#endif
