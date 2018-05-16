#include "/cosmosis/cosmosis/datablock/c_datablock.h"
#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "/cosmosis/cosmosis/datablock/section_names.h"
#include "cubacpp/cubacpp.hh"
#include "gamma_t.hh"
#include "mz_power_law.hh"

#include "test/del_sig_cen_t.hh"
#include "test/dv_do_dz_t.hh"
#include "test/ez.hh"
#include "test/ez_sq.hh"
#include "test/hmf_t.hh"
#include "test/lc_lt_t.hh"
#include "test/lc_lt_t2.hh"
#include "test/lo_lc_t.hh"
#include "test/mor_t.hh"
#include "test/primitives.hh"
#include "test/roffset_t.hh"
#include "test/zo_zt_t.hh"

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "test/interp_1d.hh"
#include "test/transform.hh"

using y3_cluster::IntegrationRange;
using y3_cluster::Interp1D;
using y3_cluster::Interp2D;
using y3_cluster::mz_power_law;

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
  explicit DEL_SIG_CEN_y1(double c) : _c(c) {}

  double
  operator()(double, double lnM, double) const
  {
    return std::exp(lnM);
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
  while (file1 >> num) {
    dndlnmh.push_back(num);
  }
  if (dndlnmh.empty())
    return 1;

  std::vector<double> mh;
  std::ifstream file2(
    "/cosmosis/cosmosis-standard-library/y3_cluster_cpp/test/m_h.txt");
  while (file2 >> num) {
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
  y3_cluster::LC_LT_t lc_lt{1.24, 4.19, 2.03, 0.32, 0.12};
  y3_cluster::ZO_ZT_t zo_zt{0.05};
  y3_cluster::ROFFSET_t roffset{0.2};
  T_CEN_t t_cen;
  T_MIS_t t_mis;
  A_CEN_t a_cen;
  A_MIS_t a_mis;
  auto p1 = std::make_shared<Interp2D const>(mh, zz, dndlnmh);
  y3_cluster::HMF_t hmf(p1, 0.037, 1.008);
  // DEL_SIG_CEN_t dsc{5., 0.5};
  y3_cluster::DEL_SIG_CEN_t dsc{5.};
  // DEL_SIG_MIS_t dsc{5., 0.5};
  DEL_SIG_MIS_t dsm;

  auto da_f = std::make_shared<Interp1D const>(zz, da_arr);
  y3_cluster::DV_DO_DZ_t dvdodz(da_f, y3_cluster::EZ(0.3, 0.7, 0));
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
