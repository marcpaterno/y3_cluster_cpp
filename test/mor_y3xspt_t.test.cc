#include "catch2/catch.hpp"

#include "models/sptxdes/mor_y3xspt_t.hh"
#include "utils/interp_1d.hh"
#include "utils/read_vector.hh"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using y3_cluster::Interp1D;
using y3_cluster::MOR_Y3xSPT_t;

TEST_CASE("MOR_Y3xSPT_t works")
{
  // Nuisance parameters
  double const lamtrue = 92.0;
  double const zeta = 10.0;
  double const ztrue = 0.40;
  double const lnm200m = 34.5;
  double const truth = 2.67235429418746436858e-03;

  // Model parameters
  double const Mmin = 1.35e11;
  double const M1 = 2.34e12;
  double const alpha = 0.748;
  double const epsilon = -0.07;
  double const sigintr = 0.30;
  double const Asz = 4.08;
  double const Bsz = 1.65;
  double const Csz = 0.64;
  double const siglnzeta = 0.20;
  double const r = 0.3;
  double const zplam = 0.45;
  double const lnMpsz = 33.33;
  double const zpsz = 0.60;
  double const omega_m = 0.30;
  double const omega_l = 0.70;
  double const omega_k = 0.00;

  // Load the mass conversion
  auto lnm200m_in = read_vector(
    "mass_conversion/lnm200m_z0.5_child18_msunhinv_mor_y3xspt_t.txt");
  auto lnm500c_in = read_vector(
    "mass_conversion/lnm500c_z0.5_child18_msunhinv_mor_y3xspt_t.txt");
  auto mconv = std::make_shared<Interp1D const>(lnm200m_in, lnm500c_in);

  // Set up the MOR object
  MOR_Y3xSPT_t mor = MOR_Y3xSPT_t(Mmin,
                                  M1,
                                  alpha,
                                  epsilon,
                                  sigintr,
                                  Asz,
                                  Bsz,
                                  Csz,
                                  siglnzeta,
                                  r,
                                  zplam,
                                  lnMpsz,
                                  zpsz,
                                  omega_m,
                                  omega_l,
                                  omega_k,
                                  mconv);

  // Test code against truth
  double constexpr epsrel = 1.0e-6;
  double constexpr epsabs = 1.0e-12;
  double const testval = mor(lamtrue, zeta, ztrue, lnm200m);
  CHECK(testval == Approx(truth).epsilon(epsrel).margin(epsabs));
}
