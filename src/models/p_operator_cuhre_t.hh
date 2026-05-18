// P_operator<WeightF> -- Costanzi-2026 P[X] operator integrand (eq. 4 of
// docs/richness_selection.tex), one CosmoSIS scalar-integration class
// template shared by three thin `.cc` drivers (P1, I1, I2).
//
// Flow mirrors richness_selection/sel_bias.py::_P_operator:
//
//   P[X](lob, zob)  =  int dlt dz dlnM   n(M, z) * dV/dOmega/dz * Omega(z)
//                                      * B_i(lt, z) * P(lt | M, z)
//                                      * w_z(z, zob)
//                                      * 2pi * int dtheta sin(theta) f_A(theta; lt, lob, zob)
//                                                                    * X(theta, ...)
//
// cubacpp::Cuhre (adaptive PAGANI-style cubature) drives the OUTER 3-D
// integral over (lt, zt, lnM).  The INNER theta integral lives inside
// operator() as a GL loop with z-dependent lower bound theta_excl(z)
// (split at exclusion).
//
// The wall grid is 1-D over lambda_bin (integer index into the 4 Y3
// richness bins), exactly like num_counts_full.  Output `vals` is
// 2-D: N_lob rows x N_zob columns (zob is the Cartesian-product axis
// via wall-of-numbers zob_low / zob_high).  Actually, per the agreement:
// wall grid = (zob_low, zob_high) x lambda_bin is cartesian-producted
// via use_cartesian_product=T, so vals.shape = (N_lob, N_zob).
//
// Template parameter WeightF selects which scalar is integrated:
//   P1: F = 1
//   I2: F = b(M,z) * xi_NL(Delta_chi)
//   I1: F = b(M,z) * xi_NL(Delta_chi) * sigmoid(theta)
#ifndef Y3_CLUSTER_CPP_P_OPERATOR_CUHRE_T_HH
#define Y3_CLUSTER_CPP_P_OPERATOR_CUHRE_T_HH

#include "utils/make_grid_points.hh"
#include "utils/make_integration_volumes.hh"

#include "cosmosis/datablock/datablock.hh"
#include "cubacpp/integration_volume.hh"

#include "utils/interp_1d.hh"
#include "utils/interp_2d.hh"
#include "utils/make_interp_1d.hh"
#include "utils/make_interp_2d.hh"
#include "models/z_kernel_data.hh"

#include "models/dv_do_dz_t.hh"
#include "models/hmf_t.hh"
#include "models/mor_hod_t.hh"
#include "models/omega_z_des.hh"
#include "models/projection_y3_b_i_t.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <vector>

namespace y3_cluster {

  namespace p_op_detail {

    inline constexpr double PI = 3.14159265358979323846;

    inline double R_lambda(double lam) { return std::pow(lam / 100.0, 0.2); }

    // GL nodes on [a, b] (Newton iteration).
    inline void
    gl_nodes(double a, double b, std::size_t N,
             std::vector<double>& xs, std::vector<double>& ws)
    {
      xs.assign(N, 0.0);
      ws.assign(N, 0.0);
      double const eps  = 1e-14;
      double const mid  = 0.5 * (a + b);
      double const hlen = 0.5 * (b - a);
      for (std::size_t i = 0; i < (N + 1) / 2; ++i) {
        double z = std::cos(PI * (double(i) + 0.75) / (double(N) + 0.5));
        double z1, pp;
        do {
          double p1 = 1.0, p2 = 0.0;
          for (std::size_t j = 0; j < N; ++j) {
            double const p3 = p2;
            p2 = p1;
            p1 = ((2.0 * double(j) + 1.0) * z * p2 - double(j) * p3) /
                 (double(j) + 1.0);
          }
          pp = double(N) * (z * p1 - p2) / (z * z - 1.0);
          z1 = z;
          z  = z1 - p1 / pp;
        } while (std::abs(z - z1) > eps);
        xs[i]         = mid - hlen * z;
        xs[N - 1 - i] = mid + hlen * z;
        double const w = 2.0 * hlen / ((1.0 - z * z) * pp * pp);
        ws[i]         = w;
        ws[N - 1 - i] = w;
      }
    }

    // Two-disk fractional overlap f_A(theta; theta_lob, theta_lt)
    // (richness_selection/geometry.py::area_overlap, closed form).
    inline double
    area_overlap(double theta, double theta_lob, double theta_lt)
    {
      if (theta > theta_lob + theta_lt) return 0.0;
      if (theta <= std::abs(theta_lob - theta_lt)) {
        // Full containment, normalised by pi * theta_lt^2 (the
        // projector area).  Matches richness_selection.geometry.area_overlap:
        //   * theta_lt < theta_lob  : projector fully inside the target, overlap = pi*theta_lt^2 -> 1
        //   * theta_lt > theta_lob  : target fully inside the projector, overlap = pi*theta_lob^2 -> (theta_lob/theta_lt)^2
        if (theta_lt <= theta_lob) return 1.0;
        double const r = theta_lob / theta_lt;
        return r * r;
      }
      double const tt = theta, ll = theta_lt, lo = theta_lob;
      double c1 = (tt*tt + ll*ll - lo*lo) / (2.0 * tt * ll);
      double c2 = (tt*tt + lo*lo - ll*ll) / (2.0 * tt * lo);
      c1 = std::clamp(c1, -1.0, 1.0);
      c2 = std::clamp(c2, -1.0, 1.0);
      double sq = (-tt + ll + lo) * (tt + ll - lo)
                * ( tt - ll + lo) * (tt + ll + lo);
      sq = std::max(sq, 0.0);
      double const A = ll*ll * std::acos(c1)
                     + lo*lo * std::acos(c2)
                     - 0.5   * std::sqrt(sq);
      return A / (PI * ll * ll);
    }

    // Bin centres of the 4 Y3 richness bins:
    //   [20, 30) [30, 45) [45, 60) [60, ~200)
    inline double
    lob_center(int bin)
    {
      // DES-Y3 richness-bin arithmetic centres for edges
      // [20, 30, 45, 60, 200]: (25, 37.5, 52.5, 130).
      static constexpr std::array<double, 4> centres{{25.0, 37.5, 52.5, 130.0}};
      return centres[bin];
    }

  }  // namespace p_op_detail


  // =====================================================================
  // Shared scalar-integrand template.  The three .cc files typedef this
  // with a WeightF that returns the theta-integrand weight.
  //
  //   WeightF::module_label()          -- unique section name per scalar
  //   WeightF::uses_bias               -- if true, multiply by b(M,z)
  //   WeightF::weight(xi_val, sigmoid) -- the theta-integrand factor
  // =====================================================================
  template <typename WeightF>
  class P_operator_cuhre {
   public:
    using grid_t       = y3_cluster::grid_t<3>;       // (zob_low, zob_high, lambda_bin)
    using grid_point_t = grid_t::value_type;

   private:
    using volume_t = cubacpp::IntegrationVolume<3>;   // (lt, zt, lnM)

    std::optional<y3_cluster::HMF_t>       hmf_;
    std::optional<Interp2D>                hmb_;   // haloModel/bias(lnM, z)
    std::optional<y3_cluster::DV_DO_DZ_t>  dv_do_dz_;
    std::optional<y3_cluster::OMEGA_Z_DES> omega_z_;
    std::optional<y3_cluster::MOR_HOD_t>   mor_;
    std::optional<Interp2D>                xi_nl_;   // direct [xi_nl] lookup
    std::array<std::optional<y3_cluster::PROJ_Y3_B_i_t>, 4> b_i_;

    // Direct datablock reads for geometry (chi(z)) and photo-z sigma(z).
    std::optional<Interp1D> chi_;       // from distances section
    std::optional<Interp1D> sigma_z_;   // from models/z_kernel_data.hh

    // Current grid point
    double zob_low_  = 0.0;
    double zob_high_ = 0.0;
    double zob_      = 0.0;
    int    current_bin_ = 0;

    // Integrand config (Nth is the inner theta GL node count)
    std::size_t Nth_;

    // HMF_t mass-shift: pass lnM + ln_mass_shift_ to recover dn/dlnM
    // at exp(lnM) in physical M_sun/h.  See set_sample.
    double ln_mass_shift_ = 0.0;
    double h0_            = 1.0;   // distances/d_c [Mpc] -> cMpc/h

   public:
    explicit P_operator_cuhre(cosmosis::DataBlock& cfg)
      : Nth_(10)
    {
      for (int i = 0; i < 4; ++i) b_i_[i].emplace(i);
      if (cfg.has_val(module_label(), "n_theta"))
        Nth_ = cfg.view<int>(module_label(), "n_theta");
    }

    void
    set_sample(cosmosis::DataBlock& sample)
    {
      hmf_.emplace(sample);
      hmb_.emplace(make_Interp2D(sample, "haloModel", "lnM", "z", "bias"));
      dv_do_dz_.emplace(sample);
      omega_z_.emplace(sample);
      mor_.emplace(sample);
      // xi_NL(r, z) -- pre-tabulated by haloModelCosmosis.py (section
      // "xi_nl"). Reading a flat 2-D table is ~1000x faster than the
      // on-the-fly QAG Hankel transform in xinl.hh.
      xi_nl_.emplace(make_Interp2D(sample, "xi_nl", "r", "z", "xi_nl"));
      // Geometry: chi(z) from distances (used for Delta_chi and theta_excl).
      chi_.emplace(make_Interp1D(sample, "distances", "z", "d_c"));
      // Photo-z sigma(z) spline (in-code data, see models/z_kernel_data.hh)
      sigma_z_.emplace(Interp1D(y3_cluster::z_kernel_z(),
                                y3_cluster::z_kernel_sigma()));
      // HMF_t x-axis is ln(m_h * (omega_m-omega_nu)), so we must pass
      // lnM + ln_mass_shift_ to recover dn/dlnM at the intended mass.
      double const omm =
          sample.view<double>("cosmological_parameters", "omega_M");
      double const omn =
          sample.view<double>("cosmological_parameters", "omega_nu");
      // HMF_t already rescales its spline axis by (omm - omn) internally
      // (hmf_t.hh::_adjust_to_log), so the shift here must stay zero.
      ln_mass_shift_ = 0.0;
      // distances/d_c is in Mpc; everything else in this module lives
      // on the cMpc/h axis.  Rescale chi(z) by h0 at query time.
      h0_ = sample.view<double>("cosmological_parameters", "h0");
    }

    void
    set_grid_point(grid_point_t const& pt)
    {
      zob_low_     = pt[0];
      zob_high_    = pt[1];
      zob_         = 0.5 * (zob_low_ + zob_high_);
      current_bin_ = static_cast<int>(pt[2]);
    }

    // Integrand at PAGANI node (lt, zt, lnM).
    double
    operator()(double lt, double zt, double lnM) const
    {
      // Photo-z kernel w_z(zt, zob): 1 - u^2 if |u|<1, else 0.
      double const s = (*sigma_z_).clamp(zt);
      double const u = (zt - zob_) / s;
      if (std::abs(u) >= 1.0) return 0.0;
      double const wz = 1.0 - u * u;

      // Outer weights that do not depend on theta.  No Omega(z) and no
      // B_i(lt, z) factors -- we follow the Python reference
      // sel_bias.py::_P_operator which takes lob as a scalar and
      // leaves the survey area out of the P[X] operator (it cancels
      // in every downstream b_sel ratio).
      double const common = (*hmf_)(lnM + ln_mass_shift_, zt) *
                            (*dv_do_dz_)(zt) * (*mor_)(lt, lnM, zt);

      // Inner theta integral
      double const lob       = p_op_detail::lob_center(current_bin_);
      double const chi_o     = (*chi_).clamp(zob_) * h0_;   // cMpc/h
      double const chi_z     = (*chi_).clamp(zt)   * h0_;   // cMpc/h
      double const theta_lob = p_op_detail::R_lambda(lob) * (1.0 + zob_) / chi_o;
      double const R_excl    = p_op_detail::R_lambda(lob) * (1.0 + zob_);

      // theta_excl(zt)
      double cos_excl = (chi_z*chi_z + chi_o*chi_o - R_excl*R_excl)
                      / (2.0 * chi_z * chi_o);
      cos_excl = std::clamp(cos_excl, -1.0, 1.0);
      double const th_lo = (cos_excl >= 1.0 - 1e-12) ? 1e-6 : std::acos(cos_excl);
      double const th_hi = 2.0 * theta_lob;
      if (th_lo >= th_hi) return 0.0;

      std::vector<double> ths, wths;
      p_op_detail::gl_nodes(th_lo, th_hi, Nth_, ths, wths);

      double const theta_lt  = p_op_detail::R_lambda(lt) * (1.0 + zt) / chi_z;
      double const k_sig     = 2.5 / theta_lob;
      double const th0_sig   = 0.5 * theta_lob;

      double th_sum = 0.0;
      for (std::size_t it = 0; it < ths.size(); ++it) {
        double const th     = ths[it];
        double const sin_th = std::sin(th);
        double const fA     = p_op_detail::area_overlap(th, theta_lob, theta_lt);
        if (fA <= 0.0) continue;

        double const dchi_3d = std::sqrt(std::max(
            chi_z*chi_z + chi_o*chi_o - 2.0 * chi_z * chi_o * std::cos(th), 0.0));
        double const xi_val  = (*xi_nl_).clamp(dchi_3d, zob_);
        double const sigmoid = 1.0 / (1.0 + std::exp(-k_sig * (th - th0_sig)));
        double const w_th    = wths[it] * 2.0 * p_op_detail::PI * sin_th;

        th_sum += w_th * fA * WeightF::weight(xi_val, sigmoid);
      }

      // rho_prj has an explicit factor of lt (the projected true richness).
      double const b_mz = WeightF::uses_bias ? hmb_->clamp(lnM, zt) : 1.0;
      return wz * lt * common * b_mz * th_sum;
    }

    static char const* module_label() { return WeightF::module_label(); }

    static std::vector<volume_t>
    make_integration_volumes(cosmosis::DataBlock& cfg)
    {
      return y3_cluster::make_integration_volumes_wall_of_numbers(
          cfg, module_label(), "lt", "zt", "lnm");
    }

    static grid_t
    make_grid_points(cosmosis::DataBlock& cfg)
    {
      return y3_cluster::make_grid_points_wall_of_numbers(
          cfg, module_label(), "zo_low", "zo_high", "lambda_bin");
    }
  };

  // ---------------------------------------------------------------------
  // Weight functors -- one per integrated scalar
  // ---------------------------------------------------------------------

  struct BSelMargP1PaganiWeight {
    static char const* module_label() { return "b_sel_marg_P1_pagani"; }
    static constexpr bool uses_bias = false;
    static double weight(double /*xi*/, double /*sigmoid*/) { return 1.0; }
  };

  struct BSelMargI2PaganiWeight {
    static char const* module_label() { return "b_sel_marg_I2_pagani"; }
    static constexpr bool uses_bias = true;
    static double weight(double xi, double /*sigmoid*/) { return xi; }
  };

  struct BSelMargI1PaganiWeight {
    static char const* module_label() { return "b_sel_marg_I1_pagani"; }
    static constexpr bool uses_bias = true;
    static double weight(double xi, double sigmoid) { return xi * sigmoid; }
  };

}  // namespace y3_cluster

#endif
