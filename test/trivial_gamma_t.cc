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
#include "test/omega_z_sdss.hh"
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

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // TEST CODE FOR DELTA SIGMA////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::vector<double> del_sig_1;
  std::ifstream file5(
    "/cosmosis/cosmosis-standard-library/y3_cluster_cpp/deltasigma/deltasigma/cosmological_parameters/deltasigma_1.txt");
  while (file5 >> num) {
    del_sig_1.push_back(num);
  }
  if (del_sig_1.empty())
    return 1;

  std::vector<double> del_sig_2;
  std::ifstream file6(
    "/cosmosis/cosmosis-standard-library/y3_cluster_cpp/deltasigma/deltasigma/cosmological_parameters/deltasigma_2.txt");
  while (file6 >> num) {
    del_sig_2.push_back(num);
  }
  if (del_sig_2.empty())
    return 1;

  std::vector<double> bm;
  std::ifstream file7(
    "/cosmosis/cosmosis-standard-library/y3_cluster_cpp/deltasigma/deltasigma/cosmological_parameters/bias.txt");
  while (file7 >> num) {
    bm.push_back(num);
  }
  if (bm.empty())
    return 1;

  std::vector<double> mh1;
  std::ifstream file8(
    "/cosmosis/cosmosis-standard-library/y3_cluster_cpp/deltasigma/deltasigma/cosmological_parameters/m_h.txt");
  while (file8 >> num) {
    mh1.push_back(std::log(num));
  }
  if (mh1.empty())
    return 1;

  std::vector<double> r_perp;
  std::ifstream file9(
    "/cosmosis/cosmosis-standard-library/y3_cluster_cpp/deltasigma/deltasigma/cosmological_parameters/r_perp.txt");
  while (file9 >> num) {
    r_perp.push_back(num);
  }
  if (r_perp.empty())
    return 1;

  std::vector<double> zz1;
  std::ifstream file10(
    "/cosmosis/cosmosis-standard-library/y3_cluster_cpp/deltasigma/deltasigma/matter_power_lin/z.txt");
  while (file10 >> num) {
    zz1.push_back(num);
  }
  if (zz1.empty())
    return 1;

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
  auto p2 = std::make_shared<Interp2D const>(r_perp, mh1, del_sig_1);
  auto p3 = std::make_shared<Interp2D const>(r_perp, zz1, del_sig_2);
  auto p4 = std::make_shared<Interp2D const>(zz1, mh1, bm);
  y3_cluster::HMF_t hmf(p1, 0.037, 1.008);
  // DEL_SIG_CEN_t dsc{5., 0.5};
  y3_cluster::DEL_SIG_CEN_t dsc(p2, p3, p4, 5.);
  // DEL_SIG_MIS_t dsc{5., 0.5};
  DEL_SIG_MIS_t dsm;

  auto da_f = std::make_shared<Interp1D const>(zz, da_arr);
  y3_cluster::DV_DO_DZ_t dvdodz(da_f, y3_cluster::EZ(0.3, 0.7, 0));
  y3_cluster::OMEGA_Z_SDSS omega_z;
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
