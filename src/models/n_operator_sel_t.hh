// N_ij operator template — number counts integrand using the precomputed
// richness-selection function S_ij(lnM, z) = S_i(lnM, z) * K_j(z) published
// by src/modules/sel_function/sel_function.py. See
// docs/richness_selection_function.tex §1 (Eq. Nij_2D):
//
//   N_ij[f] = int dlnM int dz  Omega(z) * (dV/dOmega/dz)(z) * n(M, z)
//                             * S_ij(lnM, z) * f(lnM, zt)
//
// The Sel_ij factor is a per-bin Interp2D cached once per sample; the
// integrand body is a few multiplications plus a bilinear lookup.
//
// Two specialisations:
//   NOperatorSelScalar<F>:  f(lnM, zt) -> double          (grid_t<1>)
//   NOperatorSelRadial<F>:  f(R, lnM, zt) -> double       (grid_t<2>)
//
// F contract:
//   static char const* module_label();
//   explicit F(cosmosis::DataBlock& cfg)
//   void set_sample(cosmosis::DataBlock&)
//   double operator()(lnM, zt) const    [scalar form]
//   double operator()(R, lnM, zt) const [radial form]
#ifndef Y3_CLUSTER_CPP_N_OPERATOR_SEL_T_HH
#define Y3_CLUSTER_CPP_N_OPERATOR_SEL_T_HH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "cubacpp/integration_volume.hh"

#include "utils/make_grid_points.hh"
#include "utils/make_integration_volumes.hh"

#include "models/dv_do_dz_t.hh"
#include "models/hmf_t.hh"
#include "models/omega_z_des.hh"
#include "models/sel_function_t.hh"

#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

namespace y3_cluster {

  namespace nosel_detail {
    // Read the bin count from the first extent of sel_function/S_stack.
    inline int
    n_bins_from_block(cosmosis::DataBlock& sample)
    {
      auto const& nd = sample.view<cosmosis::ndarray<double>>(
        "sel_function", "S_stack");
      return static_cast<int>(nd.extents()[0]);
    }

    // SFINAE: does the weight functor advertise a set_bin(int) hook?
    // Weights that need the current lambda-bin index (e.g. miscentering
    // weights that depend on R_lambda(lambda_ob)) opt in by defining
    //   void set_bin(int lob_bin);
    // Weights that do not advertise the hook compile to the no-op
    // overload below, so existing weight functors are undisturbed.
    template <typename W, typename = void>
    struct has_set_bin : std::false_type {};

    template <typename W>
    struct has_set_bin<
        W, std::void_t<decltype(std::declval<W&>().set_bin(0))>>
      : std::true_type {};

    template <typename W>
    inline void
    maybe_set_bin(W& w, int b, std::true_type)
    {
      w.set_bin(b);
    }

    template <typename W>
    inline void
    maybe_set_bin(W&, int, std::false_type)
    {}
  } // namespace nosel_detail


  // ======================================================================
  // Scalar-observable N operator: f(lnM, zt) -> grid_t<1> (bin_index)
  // ======================================================================
  template <typename F>
  class NOperatorSelScalar {
  public:
    using grid_t       = y3_cluster::grid_t<1>;
    using grid_point_t = grid_t::value_type;

  private:
    using volume_t = cubacpp::IntegrationVolume<2>;   // (zt, lnm)

    std::optional<y3_cluster::HMF_t>       hmf_;
    std::optional<y3_cluster::DV_DO_DZ_t>  dv_do_dz_;
    std::optional<y3_cluster::OMEGA_Z_DES> omega_z_;

    // One Interp2D cache per richness+redshift bin (built per sample by
    // sel_function.py and read here).
    std::vector<SelFunction_t> sel_;

    F weight_;
    int current_bin_{0};

  public:
    explicit NOperatorSelScalar(cosmosis::DataBlock& cfg) : weight_(cfg) {}

    void
    set_sample(cosmosis::DataBlock& sample)
    {
      hmf_.emplace(sample);
      dv_do_dz_.emplace(sample);
      omega_z_.emplace(sample);
      weight_.set_sample(sample);

      sel_.clear();
      int const n = nosel_detail::n_bins_from_block(sample);
      sel_.reserve(n);
      for (int k = 0; k < n; ++k) sel_.emplace_back(sample, k);
    }

    void
    set_grid_point(grid_point_t const& pt)
    {
      current_bin_ = static_cast<int>(pt[0]);
      nosel_detail::maybe_set_bin(weight_, current_bin_,
                                  nosel_detail::has_set_bin<F>{});
    }

    double
    operator()(double zt, double lnM) const
    {
      double const S_ij = sel_[current_bin_](lnM, zt);
      return (*hmf_)(lnM, zt) * (*dv_do_dz_)(zt) * (*omega_z_)(zt) *
             S_ij * weight_(lnM, zt);
    }

    static char const* module_label() { return F::module_label(); }

    static std::vector<volume_t>
    make_integration_volumes(cosmosis::DataBlock& cfg)
    {
      return y3_cluster::make_integration_volumes_wall_of_numbers(
        cfg, module_label(), "zt", "lnm");
    }

    static grid_t
    make_grid_points(cosmosis::DataBlock& cfg)
    {
      return y3_cluster::make_grid_points_wall_of_numbers(
        cfg, module_label(), "bin_index");
    }
  };


  // ======================================================================
  // Radial-observable N operator: f(R, lnM, zt)
  //   grid_t<2> = (bin_index, r_perp), Cartesian product
  // ======================================================================
  template <typename F>
  class NOperatorSelRadial {
  public:
    using grid_t       = y3_cluster::grid_t<2>;
    using grid_point_t = grid_t::value_type;

  private:
    using volume_t = cubacpp::IntegrationVolume<2>;   // (zt, lnm)

    std::optional<y3_cluster::HMF_t>       hmf_;
    std::optional<y3_cluster::DV_DO_DZ_t>  dv_do_dz_;
    std::optional<y3_cluster::OMEGA_Z_DES> omega_z_;

    std::vector<SelFunction_t> sel_;

    F weight_;
    int    current_bin_{0};
    double R_{0.0};

  public:
    explicit NOperatorSelRadial(cosmosis::DataBlock& cfg) : weight_(cfg) {}

    void
    set_sample(cosmosis::DataBlock& sample)
    {
      hmf_.emplace(sample);
      dv_do_dz_.emplace(sample);
      omega_z_.emplace(sample);
      weight_.set_sample(sample);

      sel_.clear();
      int const n = nosel_detail::n_bins_from_block(sample);
      sel_.reserve(n);
      for (int k = 0; k < n; ++k) sel_.emplace_back(sample, k);
    }

    void
    set_grid_point(grid_point_t const& pt)
    {
      current_bin_ = static_cast<int>(pt[0]);
      R_ = pt[1];
      nosel_detail::maybe_set_bin(weight_, current_bin_,
                                  nosel_detail::has_set_bin<F>{});
    }

    double
    operator()(double zt, double lnM) const
    {
      double const S_ij = sel_[current_bin_](lnM, zt);
      return (*hmf_)(lnM, zt) * (*dv_do_dz_)(zt) * (*omega_z_)(zt) *
             S_ij * weight_(R_, lnM, zt);
    }

    static char const* module_label() { return F::module_label(); }

    static std::vector<volume_t>
    make_integration_volumes(cosmosis::DataBlock& cfg)
    {
      return y3_cluster::make_integration_volumes_wall_of_numbers(
        cfg, module_label(), "zt", "lnm");
    }

    static grid_t
    make_grid_points(cosmosis::DataBlock& cfg)
    {
      return y3_cluster::make_grid_points_cartesian_product(
        cfg, module_label(), "bin_index", "r_perp");
    }
  };

} // namespace y3_cluster

#endif
