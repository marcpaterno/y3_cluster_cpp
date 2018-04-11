#include "cubacpp/cubacpp.hh"
#include "gamma_t.hh"
#include <chrono>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "test/interp_1d.hh"
#include "test/transform.hh"
#include <stdexcept>

using y3_cluster::Interp1D;

double constexpr pi()
{
  return 4. * std::atan(1.0);
};
double constexpr invsqrt2pi()
{
  return 1. / std::sqrt(2. * pi());
};

class mz_power_law {
public:
  mz_power_law(double A, double B, double C)
    : log2_A_(std::log2(A)), B_(B), C_(C)
  {}
  double
  operator()(double M, double z) const
  {
    double const log2_res = B_ * std::log2(M) + C_ * std::log2(1 + z) + log2_A_;
    return std::exp2(log2_res);
  }

private:
  double log2_A_;
  double B_;
  double C_;
};

inline double
gaussian(double x, double mu, double sigma)
{
  double const z = (x - mu) / sigma;
  return std::exp(-z * z / 2.) * 0.3989422804014327 / sigma;
}

class HMF_t {
public:
  HMF_t(Interp1D const* nmz, double s, double q) : _nmz(nmz), _s(s), _q(q) {}

  double
  operator()(double m, double z) const
  {
    return _nmz->eval(z) * (_s * (m - 13.8) + _q);
  }

private:
  Interp1D const* _nmz;
  double _s;
  double _q;
};

class MOR_t {
public:
  MOR_t(mz_power_law l, double sigma, double alpha)
    : _lambda(l), _sigma(sigma), _alpha(alpha)
  {}

  double
  operator()(double lt, double M, double zt) const
  {
    double const ltm = _lambda(M, zt);
    double const x = lt - ltm;
    double const erfarg = -1.0 * _alpha * (x) / (std::sqrt(2.) * _sigma);
    double const erfterm = std::erfc(erfarg);
    return gaussian(x, 0.0, _sigma) * erfterm;
  }

private:
  mz_power_law _lambda;
  double _sigma;
  double _alpha;
};

class  LO_LC_t {
public:
  explicit LO_LC_t(double alpha, double a, double b, double R_lambda) :
    _alpha(alpha), _a(a), _b(b), _R_lambda(R_lambda) {}

  double
  operator()(double lo, double lc, double R_mis) const
  {
    double x = R_mis / _R_lambda;
    double y = lo / lc;
    double mu_y = std::exp(-x*x / _alpha*_alpha);
    double sigma_y = _a * std::atan(_b*x);
    return gaussian(y, mu_y, sigma_y);
  }

private:
  double _alpha;
  double _a;
  double _b;
  double _R_lambda;
};

class LC_LT_t {
public:
  explicit LC_LT_t(double tau, double mu, double sigma, double fmsk, double fprj)
	  : _tau(tau), _mu(mu), _sigma(sigma), _fmsk(fmsk), _fprj(fprj) 
  {}

  double
  operator()(double lc, double lt, double zt) const
  {
    return (1.0-_fmsk)*(1.0-_fprj)*invsqrt2pi()*std::exp(-std::pow((lc-_mu),2.0)/(2.0*std::pow(_sigma,2.0)))/_sigma\
        +0.5*((1.0-_fmsk)*_fprj*_tau+_fmsk*_fprj/lt)*std::exp(_tau*(2.0*_mu+_tau*std::pow(_sigma,2.0)-2.0*lc)/2.0)*std::erfc((_mu+_tau*std::pow(_sigma,2.0)-lc)/(std::sqrt(2.0)*_sigma))\
        +_fmsk*(std::erfc((_mu-lc-lt)/(std::sqrt(2.0)*_sigma))-std::erfc((_mu-lc)/(std::sqrt(2.0)*_sigma)))/(2.0*lt)\
        -_fmsk*_fprj*(std::exp(-_tau*lt)*std::exp(_tau*(2.0*_mu+_tau*std::pow(_sigma,2.0)-2.0*lc)/2.0)*std::erfc((_mu+_tau*std::pow(_sigma,2.0)-lc-lt)/(std::sqrt(2.0)*_sigma)))/(2.0*lt);
  }

private:
  double _tau;
  double _mu;
  double _sigma;
  double _fmsk;
  double _fprj;
};

class ZO_ZT_t {
public:
  explicit ZO_ZT_t(double sigma) : _sigma(sigma) {}

  double
  operator()(double zo, double zt) const
  {
    double const x = zo - zt;
    return gaussian(x, 0., _sigma);
  }

private:
  double _sigma;
};

class ROFFSET_t {
public:
  explicit ROFFSET_t(double tau) : _tau(tau) {}

  double
  operator()(double x) const
  {
    // return x - offset;
    return x / _tau / _tau * std::exp(-x / _tau);
  }

private:
  double _tau;
};

struct T_CEN_t {
  double
  operator()(double x, double y) const
  {
    // return x + y;
    return 1.0;
  }
};

struct T_MIS_t {
  double
  operator()(double x, double y, double z) const
  {
    // return x + y * z;
    return 1.0;
  }
};

struct A_CEN_t {
  double
  operator()(double a, double b, double c, double d) const
  {
    // return (a + b) * (c + d);
    return 1.0;
  }
};

struct A_MIS_t {
  double
  operator()(double a, double b, double c, double d, double e) const
  {
    // return (a + b + c) * (d + e);
    return 1.0;
  }
};

struct DEL_SIG_CEN_t {
  double
  operator()(double x, double y) const
  {
    // return 3. * x + y;
    return 1.0;
  }
};

struct DEL_SIG_MIS_t {
  double
  operator()(double x, double y, double z) const
  {
    // return (2. * x + 0.5 * y) * z;
    return 1.0;
  }
};

class DV_DO_DZ_t {
public:

  double
  operator()(double x) const
  {
    return 1.0;
  }
};

class OMEGA_Z_t {
public:
  explicit OMEGA_Z_t(double zpivot) : _zpivot(zpivot) {}

  double
  operator()(double zt) const
  {
    double poly_coeff_vol[12]={ -1.14293122E05, 5.96846869E04, 9.24239180E03, -2.23118813E03, \
                                -4.52580713E03, 1.18404878E03, 1.27951911E02, -5.05716847E01, \
                                1.01744577E00,  -3.11253383E-01, 5.48481084E-03,   3.12629987E00};
    int poly_deg= 12;
    double omega_z= 0.0;

    for(int i=0; i<12; i++){
      omega_z=omega_z+ poly_coeff_vol[i]*std::pow(zt- _zpivot,poly_deg-i-1);
    }
    return omega_z;
  }

private:
  double _zpivot;
};


template <class ALG, class F>
void
time_integration(ALG alg,
                 F f,
                 double epsrel,
                 double epsabs,
                 char const* algname)
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
  std::vector<std::string> args{argv + 1, argv + argc};
  if (args.size() != 1) {
    std::cerr << "Please specify an integer maxeval\n";
  }

  std::vector<double> dndlnmh;
  double num { 0 };
  std::ifstream file1 ("/cosmosis/cosmosis-standard-library/y3_cluster_cpp/test/dndlnmh.txt");
  while (file1 >> num)
       dndlnmh.push_back(num);

  std::vector<double> mh;
  std::ifstream file2 ("/cosmosis/cosmosis-standard-library/y3_cluster_cpp/test/mh.txt");
  while (file2 >> num)
       mh.push_back(num);

  std::vector<double> zz;
  std::ifstream file3 ("/cosmosis/cosmosis-standard-library/y3_cluster_cpp/test/z.txt");
  while (file3 >> num)
       zz.push_back(num);
  if (zz.empty()) return 1;


  long long maxeval = std::stoll(args[0]);
  MOR_t mor{mz_power_law{1., 1., 0.1}, 1., 1.};
  LO_LC_t lo_lc{1.66, 0.26, 1.43, 1.0};
  LC_LT_t lc_lt{1.24, 4.19, 2.03, 0.32, 0.12};
  ZO_ZT_t zo_zt{0.1};
  ROFFSET_t roffset{2.0};
  T_CEN_t t_cen;
  T_MIS_t t_mis;
  A_CEN_t a_cen;
  A_MIS_t a_mis;
  Interp1D f{zz, zz};
  HMF_t hmf{&f, 0.037, 1.008};
  DEL_SIG_CEN_t dsc;
  DEL_SIG_MIS_t dsm;
  DV_DO_DZ_t dvdodz;
  OMEGA_Z_t omega_z;
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
                                    dsm,
				    dvdodz,
				    omega_z);
  // Call the integrand once, printing out the value.
  std::cout << gti(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9) << std::endl;
  double const epsrel = 1.0e-3;
  double const epsabs = 1.0e-12;
  cubacpp::Cuhre c;
  c.maxeval = maxeval;
  time_integration(c, gti, epsrel, epsabs, "cuhre");

  cubacpp::Vegas v;
  v.maxeval = maxeval;
  time_integration(v, gti, epsrel, epsabs, "vegas");

  cubacpp::Suave s;
  s.maxeval = maxeval;
  s.flatness = 100.0;
  time_integration(s, gti, epsrel, epsabs, "suave");
};
