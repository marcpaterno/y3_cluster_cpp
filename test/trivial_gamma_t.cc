#include "/cosmosis/cosmosis/datablock/c_datablock.h"
#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "/cosmosis/cosmosis/datablock/section_names.h"
#include "cubacpp/cubacpp.hh"
#include "gamma_t.hh"
#include "mz_power_law.hh"

#include "test/ez.hh"
#include "test/ez_sq.hh"
#include "test/hmf_t.hh"
#include "test/lo_lc_t.hh"
#include "test/mor_t.hh"
#include "test/primitives.hh"

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "test/interp_1d.hh"
#include "test/transform.hh"

using y3_cluster::IntegrationRange;
using y3_cluster::Interp1D;
using y3_cluster::mz_power_law;

class LC_LT_t2 {
public:
  double
  operator()(double lc, double lt, double /* zt */) const
  {
    return y3_cluster::gaussian(lc - lt, 0.0, 10.0);
  }
};
class LC_LT_t {
public:
  explicit LC_LT_t(double tau,
                   double mu,
                   double sigma,
                   double fmsk,
                   double fprj)
    : _tau(tau), _mu(mu), _sigma(sigma), _fmsk(fmsk), _fprj(fprj)
  {}

  explicit LC_LT_t(cosmosis::DataBlock& sample)
  {
    sample.get_val<double>("LC_LT_params", "tau", _tau);
    sample.get_val<double>("LC_LT_params", "mu", _mu);
    sample.get_val<double>("LC_LT_params", "sigma", _sigma);
    sample.get_val<double>("LC_LT_params", "fmsk", _fmsk);
    sample.get_val<double>("LC_LT_params", "fprj", _fprj);
  }

  double
  operator()(double lc, double lt, double /* zt */) const
  {
    return (1.0 - _fmsk) * (1.0 - _fprj) * y3_cluster::invsqrt2pi() *
             std::exp(-(lc - _mu) * (lc - _mu) / (2.0 * _sigma * _sigma)) /
             _sigma +
           0.5 * ((1.0 - _fmsk) * _fprj * _tau + _fmsk * _fprj / lt) *
             std::exp(_tau * (2.0 * _mu + _tau * _sigma * _sigma - 2.0 * lc) /
                      2.0) *
             std::erfc((_mu + _tau * _sigma * _sigma - lc) /
                       (std::sqrt(2.0) * _sigma)) +
           _fmsk *
             (std::erfc((_mu - lc - lt) / (std::sqrt(2.0) * _sigma)) -
              std::erfc((_mu - lc) / (std::sqrt(2.0) * _sigma))) /
             (2.0 * lt) -
           _fmsk * _fprj *
             (std::exp(-_tau * lt) *
              std::exp(_tau * (2.0 * _mu + _tau * _sigma * _sigma - 2.0 * lc) /
                       2.0) *
              std::erfc((_mu + _tau * _sigma * _sigma - lc - lt) /
                        (std::sqrt(2.0) * _sigma))) /
             (2.0 * lt);
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

  explicit ZO_ZT_t(cosmosis::DataBlock& sample)
  {
    sample.get_val<double>("ZO_ZT_params", "sigma", _sigma);
  }

  double
  operator()(double zo, double zt) const
  {
    double const x = zo - zt;
    return y3_cluster::gaussian(x, 0., _sigma);
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
  operator()(double, double) const
  {
    // return x + y;
    return 1.0;
  }
};

struct T_MIS_t {
  double
  operator()(double, double, double) const
  {
    // return x + y * z;
    return 1.0;
  }
};

struct A_CEN_t {
  double
  operator()(double, double, double, double) const
  {
    // return (a + b) * (c + d);
    return 1.0;
  }
};

struct A_MIS_t {
  double
  operator()(double, double, double, double, double) const
  {
    // return (a + b + c) * (d + e);
    return 1.0;
  }
};

/* NFW Profile , in h*M_solar*Mpc^-2 */
class DEL_SIG_CEN_y1 {
public:
  double
  operator()(double, double lnM) const
  {
    return std::exp(lnM);
  }
};

class DEL_SIG_CEN_t {
public:
  explicit DEL_SIG_CEN_t(double c) : _c(c) {}

  explicit DEL_SIG_CEN_t(cosmosis::DataBlock& sample)
  {
    sample.get_val<double>("del_sig_cen_params", "c", _c);
  }

  double
  operator()(double r, double lnM, double zt) const
  /*r in h^-1 Mpc */ /* M in h^-1 M_solar, represents M_{200} */
  {
    y3_cluster::EZ_sq ez_sq{0.3, 0.7, 0.};

    double delta_c =
      200. * _c * _c * _c / (3. * (std::log(1. + _c) - _c / (1. + _c)));
    double rho_crit = 2.77526157E11 * ez_sq(zt);
    /* EZ*EZ would be a little bit slower than direct definition */

    double r_200 = std::pow(
      3. * std::exp(lnM) / (800. * y3_cluster::pi() * rho_crit), 1. / 3.);
    double r_s = r_200 / _c;

    double r_ratio = r / r_s;
    double coeff = r_s * delta_c * rho_crit;

    if (r_ratio < 1.) {
      return coeff *
             (8. * std::atanh(std::sqrt((1. - r_ratio) / (1. + r_ratio))) /
                (r_ratio * r_ratio * std::sqrt(1. - r_ratio * r_ratio)) +
              4. * std::log(r_ratio / 2.) / (r_ratio * r_ratio) -
              2. / (r_ratio * r_ratio - 1.) +
              4. * std::atanh(std::sqrt((1. - r_ratio) / (1. + r_ratio))) /
                ((r_ratio * r_ratio - 1.) * std::sqrt(1. - r_ratio * r_ratio)));
    } else if (r_ratio == 1.) {
      return coeff * (10. / 3. + 4. * std::log(0.5));
    } else {
      return coeff *
             (8. * std::atan(std::sqrt((r_ratio - 1.) / (r_ratio + 1.))) /
                (r_ratio * r_ratio * std::sqrt(r_ratio * r_ratio - 1.)) +
              4. * std::log(r_ratio / 2.) / (r_ratio * r_ratio) -
              2. / (r_ratio * r_ratio - 1.) +
              4. * std::atan(std::sqrt((r_ratio - 1.) / (r_ratio + 1.))) /
                (std::pow(r_ratio * r_ratio - 1., 1.5)));
    }
  }

private:
  double _c;
};

struct DEL_SIG_MIS_t {
  double
  operator()(double, double, double) const
  {
    // return (2. * x + 0.5 * y) * z;
    return 1.0;
  }
};

class DV_DO_DZ_t {
public:
  DV_DO_DZ_t(Interp1D const* da, y3_cluster::EZ ezt) : _da(da), _ezt(ezt) {}

  explicit DV_DO_DZ_t(cosmosis::DataBlock& sample)
    : _da([](cosmosis::DataBlock& x) {
        std::vector<double> xs, ys;
        x.get_val("DV_D0_DZ_params", "xs", xs);
        x.get_val("DV_D0_DZ_params", "ys", ys);
        // NOTE: This is required as inputs to Interp1D must be const's
        std::vector<double> const cxs(xs), cys(ys);
	Interp1D const interp{cxs, cys};
	Interp1D const* pinterp = &interp;
	return pinterp;
      }(sample)), _ezt(EZ{1.0, 1.0, 1.0})
      //  _ezt([](cosmosis::DataBlock& x) {
      //    double omega_m, omega_l, omega_k
      //    x.get_val<double>("", "omega_m", omega_m);
      //    x.get_val<double>("", "omega_l", omega_l);
      //    x.get_val<double>("", "omega_k", omega_k);
      //    return EZ(omega_m, omega_l, omega_k);
      //}(sample))
      {}

  double
  operator()(double zt) const
  {
    double const _da_z = _da->eval(zt);
    return 3000.0 * (1.0 + zt) * (1.0 + zt) * _da_z * _da_z / _ezt(zt);
  }

private:
  Interp1D const* _da;
  y3_cluster::EZ _ezt;
};

class OMEGA_Z_t {
public:
  double
  operator()(double zt) const
  {
    double poly_coeff_vol[12] = {-1.14293122E05,
                                 5.96846869E04,
                                 9.24239180E03,
                                 -2.23118813E03,
                                 -4.52580713E03,
                                 1.18404878E03,
                                 1.27951911E02,
                                 -5.05716847E01,
                                 1.01744577E00,
                                 -3.11253383E-01,
                                 5.48481084E-03,
                                 3.12629987E00};
    int poly_deg = 12;
    double omega_z = 0.0;
    double zpivot = 0.2;

    for (int i = 0; i < 12; i++) {
      omega_z =
        omega_z + poly_coeff_vol[i] * std::pow(zt - zpivot, poly_deg - i - 1.);
    }
    return omega_z;
  }
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
  double num{0};
  std::ifstream file1(
    "/cosmosis/cosmosis-standard-library/y3_cluster_cpp/test/dndlnmh.txt");
  for (int a = 1; a < 970; a = a + 1) {
    file1 >> num;
    dndlnmh.push_back(num);
  }
  if (dndlnmh.empty())
    return 1;

  std::vector<double> mh;
  std::ifstream file2(
    "/cosmosis/cosmosis-standard-library/y3_cluster_cpp/test/m_h.txt");
  for (int a = 1; a < 970; a = a + 1) {
    file2 >> num;
    mh.push_back(std::log(num));
  }
  if (mh.empty())
    return 1;

  std::vector<double> zz;
  std::ifstream file3(
    "/cosmosis/cosmosis-standard-library/y3_cluster_cpp/test/z.txt");
  while (file3 >> num) {
    zz.push_back(num);
  }
  if (zz.empty())
    return 1;

  std::vector<double> da_arr; // in h inverse Mpc
  std::ifstream file4(
    "/cosmosis/cosmosis-standard-library/y3_cluster_cpp/test/d_a.txt");
  while (file4 >> num)
    da_arr.push_back(num);
  if (da_arr.empty())
    return 1;

  long long maxeval = std::stoll(args[0]);
  y3_cluster::MOR_t mor{mz_power_law{1.e-14, 1., 0.1}, 1., 1.};
  y3_cluster::LO_LC_t lo_lc{1.66, 0.26, 1.43, 1.0};
  LC_LT_t lc_lt{1.24, 4.19, 2.03, 0.32, 0.12};
  ZO_ZT_t zo_zt{0.05};
  ROFFSET_t roffset{0.2};
  T_CEN_t t_cen;
  T_MIS_t t_mis;
  A_CEN_t a_cen;
  A_MIS_t a_mis;
  Interp1D f{mh, dndlnmh};
  y3_cluster::HMF_t hmf{&f, 0.037, 1.008};
  // DEL_SIG_CEN_t dsc{5., 0.5};
  DEL_SIG_CEN_t dsc{5.};
  // DEL_SIG_MIS_t dsc{5., 0.5};
  DEL_SIG_MIS_t dsm;

  Interp1D da_f{zz, da_arr};
  DV_DO_DZ_t dvdodz(&da_f, y3_cluster::EZ(0.3, 0.7, 0));
  OMEGA_Z_t omega_z;
  IntegrationRange lo_ir{10, 30};
  IntegrationRange zo_ir{0.2, 0.3};
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
                                    omega_z,
                                    lo_ir,
                                    zo_ir);

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
