#include "cubacpp/cubacpp.hh"
#include "examples/gamma_t.hh"
#include <chrono>
#include <iostream>
#include <string>
#include <vector>

struct HMF_t {
  double
  operator()(double x, double y) const
  {
    return x * y;
  }
};

struct MOR_t {
  double
  operator()(double x, double y, double z) const
  {
    return x + y + z;
  }
};

struct LO_LC_t {
  double
  operator()(double x, double y, double z) const
  {
    return x * y + z;
  }
};

struct LC_LT_t {
  double
  operator()(double x, double y) const
  {
    return x + 2 * y;
  }
};

struct ZO_ZT_t {
  double
  operator()(double x, double y, double z) const
  {
    return (x + y) * z;
  }
};

class ROFFSET_t {
public:
  ROFFSET_t(double o) : offset(o) {}

  double
  operator()(double x) const
  {
    return x - offset;
  }

private:
  double offset;
};

struct T_CEN_t {
  double
  operator()(double x, double y) const
  {
    return x + y;
  }
};

struct T_MIS_t {
  double
  operator()(double x, double y, double z) const
  {
    return x + y * z;
  }
};

struct A_CEN_t {
  double
  operator()(double a, double b, double c, double d) const
  {
    return (a + b) * (c + d);
  }
};

struct A_MIS_t {
  double
  operator()(double a, double b, double c, double d, double e) const
  {
    return (a + b + c) * (d + e);
  }
};

struct DEL_SIG_CEN_t {
  double
  operator()(double x, double y) const
  {
    return 3. * x + y;
  }
};

struct DEL_SIG_MIS_t {
  double
  operator()(double x, double y, double z) const
  {
    return (2. * x + 0.5 * y) * z;
  }
};

template <class ALG, class F>
void
  time_integration(ALG alg, F f, double epsrel, double epsabs, char const* algname)
{
  auto start = std::chrono::high_resolution_clock::now();
  auto res = alg.integrate(f, epsrel, epsabs);
  auto stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = stop - start;
  std::cout << algname << ": " << res << " (" << diff.count() << "s)\n";
}

int
main(int argc, char* argv[])
{
  std::vector<std::string> args { argv+1, argv+argc };
  if (args.size() != 1) {
	std::cerr << "Please specify an integer maxeval\n";
}
  long long maxeval = std::stoll(args[0]);
  MOR_t mor;
  LO_LC_t lo_lc;
  LC_LT_t lc_lt;
  ZO_ZT_t zo_zt;
  ROFFSET_t roffset{2.0};
  T_CEN_t t_cen;
  T_MIS_t t_mis;
  A_CEN_t a_cen;
  A_MIS_t a_mis;
  HMF_t hmf;
  DEL_SIG_CEN_t dsc;
  DEL_SIG_MIS_t dsm;
  auto gti = make_gamma_t_integrand(2.0,
                                    0.11,
                                    mor,
                                    lo_lc,
                                    lc_lt,
                                    zo_zt,
                                    roffset,
                                    t_cen,
                                    t_mis,
                                    a_cen,
                                    a_mis,
                                    hmf,
                                    dsc,
                                    dsm);
  // Call the integrand once, printing out the value.
  std::cout << gti(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9) << std::endl;
  double const epsrel = 1.0e-4;
  double const epsabs = 1.0e-12;
  cubacpp::Cuhre c;
  c.maxeval = maxeval; 
  time_integration(c, gti, epsrel, epsabs, "cuhre");

  cubacpp::Vegas v;
  v.maxeval = maxeval;
  time_integration(v, gti, epsrel, epsabs, "vegas");

  cubacpp::Suave s;
  s.maxeval = maxeval;
  time_integration(s, gti, epsrel, epsabs, "suave");
};
