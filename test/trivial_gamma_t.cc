#include "/cosmosis/cosmosis/datablock/c_datablock.h"
#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "/cosmosis/cosmosis/datablock/section_names.h"
#include "cubacpp/cubacpp.hh"
#include "gamma_t.hh"
#include "mz_power_law.hh"

#include "test/a_cen_t.hh"
#include "test/a_mis_t.hh"
#include "test/del_sig_cen_t.hh"
#include "test/del_sig_mis_t.hh"
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
#include "test/t_cen_t.hh"
#include "test/t_mis_t.hh"
#include "test/zo_zt_t.hh"
#include "test/read_vector.hh"
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

// Helper template, to automate the timing of integration.
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

// The main function for exercising our integrand.
int
main(int argc, char* argv[])
{
  std::vector<std::string> args{argv + 1, argv + argc};
  if (args.size() != 1) {
    std::cerr << "Please specify an integer maxeval\n";
    return 1;
  }

  auto identity = [](double x) { return x; };
  auto log = [](double x) { return std::log(x); };

  auto const dndlnmh = read_vector("dndlnmh.txt", identity);
  auto const mh = read_vector("m_h.txt", log);
  auto const zz = read_vector("z.txt", identity);
  // da_arr in h inverse Mpc
  auto const da_arr = read_vector("d_a.txt", identity);

  auto const del_sig_1 = read_vector("deltasigma_1.txt", identity);
  auto const del_sig_2 = read_vector("deltasigma_2.txt", identity);
  auto const bm = read_vector("deltasigma_bias.txt", identity);
  auto const mh1 = read_vector("deltasigma_m_h.txt", log);
  auto const r_perp = read_vector("deltasigma_r_perp.txt", identity);
  auto const zz1 = read_vector("deltasigma_z.txt", identity);

  long long maxeval = std::stoll(args[0]);
  y3_cluster::MOR_t mor{mz_power_law{1.e-14, 1., 0.1}, 1., 1.};
  y3_cluster::LO_LC_t lo_lc{1.66, 0.26, 1.43, 1.0};
  y3_cluster::LC_LT_t lc_lt;
  y3_cluster::ZO_ZT_t zo_zt{0.05};
  y3_cluster::ROFFSET_t roffset{0.2};
  y3_cluster::T_CEN_t t_cen;
  y3_cluster::T_MIS_t t_mis;
  y3_cluster::A_CEN_t a_cen;
  y3_cluster::A_MIS_t a_mis;
  auto p1 = std::make_shared<Interp2D const>(mh, zz, dndlnmh);
  auto p2 = std::make_shared<Interp2D const>(r_perp, mh1, del_sig_1);
  auto p3 = std::make_shared<Interp2D const>(r_perp, zz1, del_sig_2);
  auto p4 = std::make_shared<Interp2D const>(zz1, mh1, bm);
  y3_cluster::HMF_t hmf(p1, 0.037, 1.008);
  // DEL_SIG_CEN_t dsc{5., 0.5};
  y3_cluster::DEL_SIG_CEN_t dsc(p2, p3, p4/*,5*/);
  // DEL_SIG_MIS_t dsc{5., 0.5};
  /*y3_cluster::DEL_SIG_MIS_t dsm;*/

  auto da_f = std::make_shared<Interp1D const>(zz, da_arr);
  y3_cluster::DV_DO_DZ_t dvdodz(da_f, y3_cluster::EZ(0.3, 0.7, 0), 0.7);
  y3_cluster::OMEGA_Z_SDSS omega_z;
  IntegrationRange lo_ir{10, 30};
  IntegrationRange zo_ir{0.2, 0.3};
  auto gti = make_gamma_t_integrand(0.7,
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
                                    dvdodz,
                                    omega_z,
                                    lo_ir,
                                    zo_ir);

  double const epsrel = 1.0e-3;
  double const epsabs = 1.0e-12;

  cubacpp::Cuhre cc;
  cc.maxeval = maxeval;
  // Won't allow integrating gti.centered directly :(
  time_integration(cc,
                   [&gti](double a, double b, double c,
                          double d, double e, double f) {
                        return gti.centered(a, b, c, d, e, f);
                   },
                   epsrel, epsabs, "cuhre");

  cubacpp::Cuhre cm;
  cm.maxeval = maxeval;
  // same deal as above
  time_integration(cm,
                   [&gti](double a, double b, double c,
                          double d, double e, double f,
                          double g, double h, double i) {
                        return gti.miscentered(a, b, c, d, e, f, g, h, i);
                   },
                   epsrel, epsabs, "cuhre");

  /*
  cubacpp::Vegas v;
  v.maxeval = maxeval;
  time_integration(v, gti, epsrel, epsabs, "vegas");

  cubacpp::Suave s;
  s.maxeval = maxeval;
  s.flatness = 100.0;
  time_integration(s, gti, epsrel, epsabs, "suave");
  */
};
