#include "/cosmosis/cosmosis/datablock/c_datablock.h"
#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "/cosmosis/cosmosis/datablock/section_names.h"
#include "cubacpp/cubacpp.hh"
#include "gamma_t.hh"
#include "mz_power_law.hh"

/*
#include "test/a_cen_t.hh"
#include "test/a_mis_t.hh"
#include "test/del_sig_y1.hh"
#include "test/del_sig_t.hh"
#include "test/dv_do_dz_t.hh"
#include "test/ez.hh"
#include "test/ez_sq.hh"
#include "test/hmf_t.hh"
#include "test/lc_lt_t.hh"
#include "test/lo_lc_t.hh"
#include "test/mor_t2.hh"
#include "test/omega_z_sdss.hh"
#include "test/primitives.hh"
#include "test/roffset_t.hh"
#include "test/t_cen_t.hh"
#include "test/t_mis_t.hh"
#include "test/zo_zt_t.hh"
#include "test/read_vector.hh"
#include "test/default_models.hh"
*/

#include <chrono>
#include <cmath>
#include <fstream>
#include <gperftools/profiler.h>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <default_models.hh>
#include <mor_t2.hh>
#include <interp_1d.hh>
#include <interp_2d.hh>
#include <read_vector.hh>
#include <transform.hh>

using y3_cluster::IntegrationRange;
using y3_cluster::Interp1D;
using y3_cluster::Interp2D;
using y3_cluster::make_gamma_t_integrand;
using y3_cluster::mz_power_law;

// Helper template, to automate the timing of integration.
template <class F>
auto
time_integration(F f,
                 char const* algname) -> decltype(f())
{
  auto start = std::chrono::high_resolution_clock::now();
  auto res = f();
  auto stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = stop - start;
  std::cout << algname << ": " << res << " (" << diff.count() << "s)\n";
  return res;
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
  // Make sure CUBA does not fork processes!
  cubacores(0, 0);

  // ============ Cosmological Parameters ============
  //            (to be passed by CosmoSIS)
  const double Omega_M = 1.87518978e-01;
  const double Omega_L = 1.0 - Omega_M;
  const double Omega_K = 0.0;
  const double h = 7.71358152e-01;

  // ============ Scaling Functions ============
  auto identity = [](double x) { return x; };
  auto log = [](double x) { return std::log(x); };
  auto log_omega_m = [Omega_M](double x) { return std::log(x*Omega_M); };

  // ============ Input Data Tables ============
  // dndlnmh.txt, m_h.txt, z.txt came from the cosmosis tinker_mf_module.so
  // d_a.txt, z_da.txt came from the cosmosis camb.so
  auto const dndlnmh = read_vector("dndlnmh.txt", identity);
  // m_h.txt is in units of: 
  //    \Omega_M M_{solar} h^{-1}
  // So, need to divide by \Omega_M to get M_{solar} h^{-1} values.
  // NOTE: 0.3 was the \Omega_M used to generate the tables, so different \Omega_M values would require different tables
  // dndlnmh is in unit of (h^3 Mpc^{-3})
  auto mh = read_vector("m_h.txt", log_omega_m);
  auto const zz = read_vector("z.txt", identity);
  // da_arr in Mpc
  auto const zz_da = read_vector("z.txt", identity);
  auto const da_arr = read_vector("d_a.txt", identity);

  auto const del_sig_1 = read_vector("deltasigma_1.txt", identity);
  auto const del_sig_2 = read_vector("deltasigma_2.txt", identity);
  auto const bm = read_vector("deltasigma_bias.txt", identity);
  auto const mh1 = read_vector("deltasigma_m_h.txt", log);
  auto const r_perp = read_vector("deltasigma_r_perp.txt", identity);
  auto const zz1 = read_vector("deltasigma_z.txt", identity);

  // ============ Integral Components ============
  // Create each term which will comprise the gamma_t integral
  // TODO: remove magic numbers
  struct MODELS : public y3_cluster::DefaultModels {
      using MOR = y3_cluster::MOR_t2;
  };
  long long maxeval = std::stoll(args[0]);
  double sigma_intr=1.29339555e-01 ;//this is a parameter that should come from cosmosis
  double alpha=6.91589257e-01 ;//this is a parameter that should come from cosmosis
  MODELS::MOR mor{pow(10,1.11375214e+01), pow(10,12.4225835912), alpha, sigma_intr};
  MODELS::LO_LC lo_lc{1.66, 0.26, 1.43, 1.0};
  MODELS::LC_LT lc_lt;
  MODELS::ZO_ZT zo_zt{0.005};
  MODELS::ROFFSET roffset{0.2};
  MODELS::T_CEN t_cen;
  MODELS::T_MIS t_mis;
  MODELS::A_CEN a_cen;
  MODELS::A_MIS a_mis;

  auto p1 = std::make_shared<Interp2D const>(mh, zz, dndlnmh);
  auto p2 = std::make_shared<Interp2D const>(r_perp, mh1, del_sig_1);
  auto p3 = std::make_shared<Interp2D const>(r_perp, zz1, del_sig_2);
  auto p4 = std::make_shared<Interp2D const>(zz1, mh1, bm);
  auto da_f = std::make_shared<Interp1D const>(zz_da, da_arr);

  MODELS::HMB bmz;
  MODELS::HMF hmf(p1, 4.50732047e-02, 1.01958078e+00);
  // TODO: Change to DEL_SIG_Y1
  // MODELS::DEL_SIG ds(p2, p3, p4);
  MODELS::DEL_SIG ds;
  // dvdodz in unit of h^{-3} Mpc^3, note that da_arr needs to be in unit of Mpc
  MODELS::DV_DO_DZ dvdodz(da_f, y3_cluster::EZ(Omega_M, Omega_L, Omega_K), h);
  MODELS::OMEGA_Z omega_z;

  IntegrationRange lo_ir{20, 27.9};
  IntegrationRange zo_ir{0.1, 0.3};

  auto gti = make_gamma_t_integrand<MODELS, 10>(0.7,
                                    mor,
                                    lo_lc,
                                    lc_lt,
                                    zo_zt,
                                    roffset,
                                    t_cen,
                                    t_mis,
                                    a_cen,
                                    a_mis,
                                    bmz,
                                    hmf,
                                    ds,
                                    dvdodz,
                                    omega_z,
                                    {lo_ir},//, {27.9, 37.6}, {37.6, 50.3}, {50.3, 69.3}},
                                    {zo_ir});

  // ============ Actual Integrations ============
  // Integrate centered and miscentered, simultaneously over all bins,
  // then each bin individually
  double const epsrel = 1.0e-3;
  double const epsabs = 1.0e-12;

  ProfilerStart("/cosmosis/cosmosis-standard-library/y3_cluster_cpp/dump.txt");
  cubacpp::Cuhre c;
  c.maxeval = maxeval;
  time_integration([&]() { return gti.integrate_centered(c, epsrel, epsabs).first; },
                   "centered-cuhre");

  time_integration([&]() { return gti.integrate_miscentered(c, epsrel, epsabs).first; },
                   "miscentered-cuhre");

  ProfilerFlush();
  ProfilerStop();
};
