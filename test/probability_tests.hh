#include "catch2/catch.hpp"
#include "cubacpp/cubacpp.hh"
#include "models/lc_lt_t.hh"
#include "models/lo_lc_t.hh"
#include "models/mor_sdss_t.hh"
#include "models/mor_t.hh"
#include "models/roffset_t.hh"
#include "utils/integration_range.hh"

#include "param_space_explorer.hh"
#include <iostream>

/******************************************
 *      Probability Integration Tests
 ******************************************
 *
 * For any dependent probability P(x|A,B,...), the following should hold:
 *
 *      \int_{xrange} P(x|A,B,...) dx = 1
 *
 * The total probability should, naturally be 1.
 *
 *   The following functions test that this holds for each probability
 * distribution, across a range of A, B, etc. values.
 *
 *   Naturally, the above equality cannot be expected to hold _exactly_, due to
 * both the precision of the integrator and the fact that the `xrange`s used
 * are generally not the complete range of possible x values. On account of
 * this each function uses a different acceptable margin.
 *
 *   Each of the below functions take the boolean flags `print` and `test`. If
 * `test` is true, it will run catch2 tests, if `print` is true, it will output
 * a CSV formatted table of results of each integration.
 */

/* TODO:
 * These argument lists are starting to get byzantine - may be simpler to pass
 * a Cosmosis data block with everything specified instead.
 */

namespace y3_cluster {

  /* Computes:
   *      \int_{lc_ir} P(\lc|\lt, \zt) d\lc
   *
   * Using `lt_width` distinct \lt values, and `zt_width` distinct \zt values in
   * the range `lt_ir` and `zt_ir`, respectively.
   */
  template <typename Integrator>
  void
  test_integrate_lc_lt(Integrator I,
                       IntegrationRange lc_ir,
                       IntegrationRange lt_ir,
                       IntegrationRange zt_ir,
                       std::size_t lt_width,
                       std::size_t zt_width,
                       const double epsrel = 1.0e-4,
                       const double epsabs = 1.0e-12,
                       bool print = false,
                       bool test = true)
  {
    if (print) std::cout << "lt,zt,lc_lt_integrated,status,error,prob\n";

    LC_LT_t lc_lt;
    double lt, zt;
    // Vary lt and zt, in the ranges lt_ir and zt_ir, with lt_width and
    // zt_width distinct values of each
    ParamSpaceExplorer<2> pse{{&lt, &zt}, {lt_ir, zt_ir}, {lt_width, zt_width}};

    do {
      const auto res = I.integrate(
        [lt, zt, lc_ir, lc_lt](double scaled_lc) {
          const double lc = lc_ir.transform(scaled_lc);
          return lc_ir.jacobian() * lc_lt(lc, lt, zt);
        },
        epsrel,
        epsabs);

      if (print)
        std::cout << lt << ", " << zt << ", " << res.value << ", " << res.status
                  << ", " << res.error << ", " << res.prob << '\n';

      if (test) {
        // In reality - since we are integrating [1, 200] not [0, +inf] - it
        // won't be exactly 1.0. So, arbitrary wiggle room 1e-2
        // TODO: should this epsilon be changed?
        CHECK(res.status == 0);
        CHECK(res.value == Approx(1.0).epsilon(1e-2).margin(1e-2));
      }
    } while (pse.next());
  }

  /* Computes:
   *      \int_{lo_ir} P(\lo|\lc, R) d\lo
   *
   * Using `lc_width` distinct \lc values, and `R_width` distinct R values in
   * the range `lc_ir` and `R_ir`, respectively.
   */
  template <typename Integrator>
  void
  test_integrate_lo_lc(Integrator I,
                       IntegrationRange lo_ir,
                       IntegrationRange lc_ir,
                       IntegrationRange R_ir,
                       std::size_t lc_width,
                       std::size_t R_width,
                       const double epsrel = 1.0e-4,
                       const double epsabs = 1.0e-12,
                       bool print = false,
                       bool test = true)
  {
    if (print) std::cout << "lc,R,lc_lt_integrated,status,error,prob\n";

    LO_LC_t lo_lc{1.66, 0.26, 1.43, 1.0};

    double lc, R;
    ParamSpaceExplorer<2> pse{{&lc, &R}, {lc_ir, R_ir}, {lc_width, R_width}};

    do {
      const auto res = I.integrate(
        [lc, R, lo_ir, lo_lc](double scaled_lo) {
          const double lo = lo_ir.transform(scaled_lo);
          return lo_ir.jacobian() * lo_lc(lo, lc, R);
        },
        epsrel,
        epsabs);

      if (print)
        std::cout << lc << ", " << R << ", " << res.value << ", " << res.status
                  << ", " << res.error << ", " << res.prob << '\n';

      if (test) {
        CHECK(res.status == 0);
        CHECK(res.value == Approx(1.0).epsilon(1e-2).margin(5e-2));
      }
    } while (pse.next());
  }

  /* Computes:
   *      \int_{R_ir} P(R_{mis}) dR_{mis}
   *
   * The only parameter is the width value `tau`. Only one `tau` value is used.
   */
  template <typename Integrator>
  void
  test_integrate_roffset(Integrator I,
                         IntegrationRange R_ir,
                         double tau,
                         const double epsrel = 1.0e-4,
                         const double epsabs = 1.0e-12,
                         bool print = false,
                         bool test = true)
  {
    if (print) std::cout << "roffset_integrand,status,error,prob\n";

    ROFFSET_t roffset(tau);

    const auto res = I.integrate(
      [R_ir, roffset](double scaled_R) {
        const double R = R_ir.transform(scaled_R);
        return R_ir.jacobian() * roffset(R);
      },
      epsrel,
      epsabs);

    if (print)
      std::cout << res.value << ", " << res.status << ", " << res.error << ", "
                << res.prob << '\n';

    if (test) {
      CHECK(res.status == 0);
      CHECK(res.value == Approx(1.0).epsilon(1e-2).margin(5e-2));
    }
  }

  /* Computes:
   *      \int_{lt_ir} P(\lt| ln(M), zt) d\lt
   *
   * Using `lnM_width` distinct ln(M) values, and `zt_width` distinct zt values
   * in the range `lnM_ir` and `zt_ir`, respectively.
   */
  template <typename Integrator>
  void
  test_integrate_mor(Integrator I,
                     IntegrationRange lt_ir,
                     IntegrationRange lnM_ir,
                     IntegrationRange zt_ir,
                     std::size_t lnM_width,
                     std::size_t zt_width,
                     const double epsrel = 1.0e-4,
                     const double epsabs = 1.0e-12,
                     bool print = false,
                     bool test = true)
  {
    if (print) std::cout << "lnM,M,zt,mor_integrated,status,error,prob\n";

    MOR_t mor{mz_power_law{9.1e-9, 0.65, 0.1}, 0.15, 0.65};

    double lnM, zt;
    ParamSpaceExplorer<2> pse{
      {&lnM, &zt}, {lnM_ir, zt_ir}, {lnM_width, zt_width}};
    do {
      const auto res = I.integrate(
        [lnM, zt, lt_ir, mor](double scaled_lt) {
          const double lt = lt_ir.transform(scaled_lt);
          return lt_ir.jacobian() * mor(lt, lnM, zt);
        },
        epsrel,
        epsabs);

      if (print)
        std::cout << lnM << ", " << std::exp(lnM) << ", " << zt << ", "
                  << res.value << ", " << res.status << ", " << res.error
                  << ", " << res.prob << '\n';

      if (test) {
        CHECK(res.status == 0);
        CHECK(res.value == Approx(1.0).epsilon(1e-2).margin(0.15));
      }
    } while (pse.next());
  }

  template <typename Integrator>
  void
  test_integrate_mor2(Integrator I,
                      IntegrationRange lt_ir,
                      IntegrationRange lnM_ir,
                      IntegrationRange zt_ir,
                      std::size_t lnM_width,
                      std::size_t zt_width,
                      const double epsrel = 1.0e-4,
                      const double epsabs = 1.0e-12,
                      bool print = false,
                      bool test = true)

  {
    if (print) std::cout << "lnM,M,zt,mor2_integrated,status,error,prob\n";

    MOR_sdss mor2{pow(10, 1.11375214e+01), pow(10, 12.4225835912), 0.65, 0.15};

    for (auto i = 0u; i < lnM_width; i++) {
      for (auto j = 0u; j < zt_width; j++) {
        const double lnM = lnM_ir.transform((i + 1) / ((double)lnM_width + 1));
        const double zt = zt_ir.transform((j + 1) / ((double)zt_width + 1));

        const auto res = I.integrate(
          [lnM, zt, lt_ir, mor2](double scaled_lt) {
            const double lt = lt_ir.transform(scaled_lt);
            return lt_ir.jacobian() * mor2(lt, lnM, zt);
          },
          epsrel,
          epsabs);

        if (print)
          std::cout << lnM << ", " << std::exp(lnM) << ", " << zt << ", "
                    << res.value << ", " << res.status << ", " << res.error
                    << ", " << res.prob << '\n';

        if (test) {
          CHECK(res.status == 0);
          CHECK(res.value == Approx(1.0).epsilon(1e-2).margin(0.01));
        }
      }
    }
  }
}
