#include "/cosmosis/cosmosis/datablock/c_datablock.h"
#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include "/cosmosis/cosmosis/datablock/section_names.h"
#include "cubacpp/cubacpp.hh"
#include "gamma_t.hh"
#include "mz_power_law.hh"

#include "test/a_cen_t.hh"
#include "test/a_mis_t.hh"
#include "test/del_sig_cen_t.hh"
#include "test/del_sig_cen_y1.hh"
#include "test/del_sig_mis_t.hh"
#include "test/dv_do_dz_t.hh"
#include "test/ez.hh"
#include "test/ez_sq.hh"
#include "test/hmf_t.hh"
#include "test/lc_lt_t.hh"
#include "test/lc_lt_t2.hh"
#include "test/lo_lc_t.hh"
#include "test/mor_t2.hh"
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
auto
time_integration(ALG alg,
                 F f,
                 double epsrel,
                 double epsabs,
                 char const* algname) -> decltype(alg.integrate(f, epsrel, epsabs))
{
  auto start = std::chrono::high_resolution_clock::now();
  auto res = alg.integrate(f, epsrel, epsabs);
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
  long long maxeval = std::stoll(args[0]);
  double sigma_intr=1.29339555e-01 ;//this is a parameter that should come from cosmosis
  double alpha=6.91589257e-01 ;//this is a parameter that should come from cosmosis
  //y3_cluster::MOR_t mor{mz_power_law{8.8e-9, alpha, 0.0}, sigma_intr, alpha};
  y3_cluster::MOR_t2 mor{pow(10,1.11375214e+01), pow(10,12.4225835912), alpha, sigma_intr};
  y3_cluster::LO_LC_t lo_lc{1.66, 0.26, 1.43, 1.0};
  y3_cluster::LC_LT_t lc_lt;
  y3_cluster::ZO_ZT_t zo_zt{0.005};
  y3_cluster::ROFFSET_t roffset{0.2};
  y3_cluster::T_CEN_t t_cen;
  y3_cluster::T_MIS_t t_mis;
  y3_cluster::A_CEN_t a_cen;
  y3_cluster::A_MIS_t a_mis;
  auto p1 = std::make_shared<Interp2D const>(mh, zz, dndlnmh);
  auto p2 = std::make_shared<Interp2D const>(r_perp, mh1, del_sig_1);
  auto p3 = std::make_shared<Interp2D const>(r_perp, zz1, del_sig_2);
  auto p4 = std::make_shared<Interp2D const>(zz1, mh1, bm);
  y3_cluster::HMF_t hmf(p1, 4.50732047e-02, 1.01958078e+00);
  //y3_cluster::DEL_SIG_CEN_t dsc(p2, p3, p4);
  y3_cluster::DEL_SIG_CEN_y1 dsc; // this is using y1 observable

  auto da_f = std::make_shared<Interp1D const>(zz_da, da_arr);
  y3_cluster::DV_DO_DZ_t dvdodz(da_f, y3_cluster::EZ(Omega_M, Omega_L, Omega_K), h); 
  // dvdodz in unit of h^{-3} Mpc^3, note that da_arr needs to be in unit of Mpc 
  y3_cluster::OMEGA_Z_SDSS omega_z;
  IntegrationRange lo_ir1{20, 27.9};
  IntegrationRange zo_ir1{0.1, 0.3};
  using MODELS = Models<decltype(mor),
                        decltype(lo_lc),
                        decltype(lc_lt),
                        decltype(zo_zt),
                        decltype(roffset),
                        decltype(t_cen),
                        decltype(t_mis),
                        decltype(a_cen),
                        decltype(a_mis),
                        decltype(hmf),
                        decltype(dsc),
                        decltype(dvdodz),
                        decltype(omega_z)>;
  auto gti = make_gamma_t_integrand<MODELS, 10, 4, 1>(0.7,
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
                                    {lo_ir1, {27.9, 37.6}, {37.6, 50.3}, {50.3, 69.3}},
                                    {zo_ir1});

  // ============ Actual Integrations ============
  // Integrate centered and miscentered, simultaneously over all bins,
  // then each bin individually
  double const epsrel = 1.0e-3;
  double const epsabs = 1.0e-12;

  cubacpp::Cuhre c;
  c.maxeval = maxeval;
  // Won't allow integrating gti.centered directly :(
  auto start = std::chrono::high_resolution_clock::now();
  auto centered_res = c.integrate([&gti](double a, double b, double c,
                                         double d, double e) {
                                      return gti.centered(a, b, c, d, e);
                                 },
                                 epsrel, epsabs);
  auto stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> diff = stop - start;

  auto centered_bins = make_gamma_t_integrated_bins(gti, centered_res);

  std::cout << "****** CENTERED ******\n"
            << "** status: " << centered_res.status << '\n'
            << "** simultaneous time: " << diff.count() << '\n';
  double time_count = 0.0;
  for (auto& bin : centered_bins) {
      std::array<y3_cluster::IntegrationRange, 1> new_lir = {bin.lo_ir};
      std::array<y3_cluster::IntegrationRange, 1> new_zir = {bin.zo_ir};
      auto gti_single = gti.with_bins(new_lir, new_zir);

      start = std::chrono::high_resolution_clock::now();
      auto single_res = c.integrate([&gti_single](double a, double b, double c,
                                                    double d, double e) {
                                                 return gti_single.centered(a, b, c, d, e);
                                            },
                                            epsrel, epsabs);
      stop = std::chrono::high_resolution_clock::now();
      diff = stop - start;

      auto single_bin = make_gamma_t_integrated_bins(gti_single, single_res)[0];

      std::cout << "(single - bin) status = " << single_res.status
                << ", time = " << diff.count() << '\n';

      std::cout << "(simultaneous) lo = " << bin.lo_ir << ", zo = " << bin.zo_ir
                << ", (N, Nw) = (" << bin.N << " +/- " << bin.N_error << ", "
                << bin.Nw << " +/- " << bin.Nw_error << ")\n";

      std::cout << "(single - bin) lo = " << single_bin.lo_ir << ", zo = " << single_bin.zo_ir
                << ", (N, Nw) = (" << single_bin.N << " +/- " << single_bin.N_error << ", "
                << single_bin.Nw << " +/- " << single_bin.Nw_error << ")\n";

      time_count += diff.count();
  }
  std::cout << "** total single time: " << time_count << '\n';

  // same deal as above
  start = std::chrono::high_resolution_clock::now();
  auto miscentered_res = c.integrate([&gti](double a, double b, double c,
                                            double d, double e, double f,
                                            double g, double h) {
                                         return gti.miscentered(a, b, c, d, e, f, g, h);
                                    },
                                    epsrel, epsabs);
  stop = std::chrono::high_resolution_clock::now();
  diff = stop - start;

  auto miscentered_bins = make_gamma_t_integrated_bins(gti, miscentered_res);

  std::cout << "****** MISCENTERED ******\n"
            << "** status: " << centered_res.status << '\n'
            << "** simultaneous time: " << diff.count() << '\n';
  time_count = 0.0;
  for (auto& bin : miscentered_bins) {
      std::array<y3_cluster::IntegrationRange, 1> new_lir = {bin.lo_ir};
      std::array<y3_cluster::IntegrationRange, 1> new_zir = {bin.zo_ir};
      auto gti_single = gti.with_bins(new_lir, new_zir);

      start = std::chrono::high_resolution_clock::now();
      auto single_res = c.integrate([&gti_single](double a, double b, double c,
                                                    double d, double e, double f,
                                                    double g, double h) {
                                                 return gti_single.miscentered(a, b, c, d, e, f, g, h);
                                            },
                                            epsrel, epsabs);
      stop = std::chrono::high_resolution_clock::now();
      diff = stop - start;

      auto single_bin = make_gamma_t_integrated_bins(gti_single, single_res)[0];

      std::cout << "(single - bin) status = " << single_res.status
                << ", time = " << diff.count() << '\n';

      std::cout << "(single - bin) lo = " << single_bin.lo_ir << ", zo = " << single_bin.zo_ir
                << ", (N, Nw) = (" << single_bin.N << " +/- " << single_bin.N_error << ", "
                << single_bin.Nw << " +/- " << single_bin.Nw_error << ")\n";

      std::cout << "(simultaneous) lo = " << bin.lo_ir << ", zo = " << bin.zo_ir
                << ", (N, Nw) = (" << bin.N << " +/- " << bin.N_error << ", "
                << bin.Nw << " +/- " << bin.Nw_error << ")\n";

      time_count += diff.count();
  }
  std::cout << "** total single time: " << time_count << '\n';
};
