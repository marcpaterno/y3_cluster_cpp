#ifndef Y3_CLUSTER_GAMMA_T_HH
#define Y3_CLUSTER_GAMMA_T_HH

#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include <sigma_crit_inverse_t.hh>
#include <integration_range.hh>
#include <read_vector.hh>
#include <transform.hh>

#include <algorithm>
#include <array>
#include <cubacpp/cubacpp.hh>
#include <cmath>
#include <iostream>
#include <string>

// The most up-to-date notes on this code are here:
// on https://v1.overleaf.com/18758665bfvkhrtbqnvz#/70551749/.

namespace y3_cluster {

// Forward declaration of our integrand class
template <typename MODELS>
class Gamma_T_Integrand;

/* The integration results of a particular (richness, redshift) bin.
 *
 * The bin is specified by `lo_ir` and `zo_ir` for richness, redshift, respectively.
 * Integrated values are in `gamma_ts`, `N`, `Nw`, `Nb`, with error and probability
 * values for each in `*_error[s]` and `*_prob[s]`.
 */
struct Gamma_T_Integrated_Bin_Result {
  y3_cluster::IntegrationRange lo_ir, zo_ir;
  std::vector<double>  radius;

  // Separate gamma_t values/errors/probs for each source bin
  std::vector<std::vector<double>>  gamma_ts;
  std::vector<std::vector<double>>  gamma_t_errors;
  std::vector<std::vector<double>>  gamma_t_probs;

  double N, N_error, N_prob,
         Nw, Nw_error, Nw_prob,
         Nb, Nb_error, Nb_prob;

  Gamma_T_Integrated_Bin_Result() : lo_ir{0.0, 1.0}, zo_ir{0.0, 1.0} {}

  template<typename MODELS>
  Gamma_T_Integrated_Bin_Result(std::size_t which_richness,
                                std::size_t which_redshift,
                                const Gamma_T_Integrand<MODELS> &gt,
                                const cubacpp::integration_results_v &results)
      : lo_ir(gt.lo_ir_[which_richness])
      , zo_ir(gt.zo_ir_[which_redshift])
      , radius(gt.r)
      , gamma_ts()
      , gamma_t_errors()
      , gamma_t_probs()
  {
    auto const nradii = radius.size(),
               npzsource = gt.npzsource;

    // The location of these bin results, as an offset in the integration output
    const auto base = (which_richness * gt.zo_ir_.size () + which_redshift)
                    * (nradii * npzsource + 3);

    // First: obtain the various number count results
    const auto nbase = base + nradii*npzsource;
    N = results.value[nbase];
    Nw = results.value[nbase + 1];
    Nb = results.value[nbase + 2];

    N_error = results.error[nbase];
    Nw_error = results.error[nbase + 1];
    Nb_error = results.error[nbase + 2];

    N_prob = results.prob[nbase];
    Nw_prob = results.prob[nbase + 1];
    Nb_prob = results.prob[nbase + 2];

    // For each WL source bin: collect gamma_t values
    for (auto i = 0u; i < npzsource; i++) {
      std::vector<double> gamma_t_src_bin;
      std::vector<double> errors_src_bin;
      std::vector<double> probs_src_bin;

      // At each radius from the cluster center: get the <gamma_t> value
      auto const zsrc_base = base + i * nradii;
      for (auto j = 0u; j < nradii; j++) {
        // Integration result is a "weighted sum" over all clusters, must
        // divide by number of clusters to obtain absolute gamma_t.
        // If we expected to see ~0 clusters, avoid divide-by-zero by making
        // gamma_t zero.
        if (Nw < 1e-2) {
          gamma_t_src_bin.push_back(0.0);
          errors_src_bin.push_back(0.0);
        } else {
          gamma_t_src_bin.push_back(results.value[zsrc_base + j] / Nw);
          errors_src_bin.push_back(results.error[zsrc_base + j] / Nw);
        }
        probs_src_bin.push_back(results.prob[zsrc_base + j]);
      }

      gamma_ts.push_back(gamma_t_src_bin);
      gamma_t_errors.push_back(errors_src_bin);
      gamma_t_probs.push_back(probs_src_bin);
    }
  }
};


/* The integration results for a collection of (richness, redshift) bins.
 *
 * The integration bins all share the `neval`, `nregions`, and `status` parameters
 * from the integration algorithm (see `cubacpp` for details).
 *
 * Results for each specific bin can be accessed from `Gamma_T_Integrated_Bin_Result_S`
 * as in a vector.
 */
struct Gamma_T_Integrated_Bin_Result_S
  : std::vector <Gamma_T_Integrated_Bin_Result>
{
  long long neval;
  int nregions = -1;
  int status = 1;

  std::size_t const n_richness;
  std::size_t const n_redshift;

  Gamma_T_Integrated_Bin_Result_S() = delete;

  Gamma_T_Integrated_Bin_Result_S (std::size_t i, std::size_t e, long long neval, int nregions, int status)
    : std::vector<Gamma_T_Integrated_Bin_Result>  (i*e)
    , neval {neval}
    , nregions {nregions}
    , status {status}
    , n_richness {i}
    , n_redshift {e}
  {}
};

inline std::ostream&
operator<<(std::ostream& os, Gamma_T_Integrated_Bin_Result_S const& res)
{
  os << "neval: " << res.neval << " nregions: " << res.nregions
     << " status: " << res.status << '\n';
  for (auto const& bin : res) {
    // Print out bin info
    os << "Bin: [zmin, zmax] = ["
       << bin.zo_ir.transform(0.0) << ", " << bin.zo_ir.transform(1.0)
       << "], [lomin, lomax] = ["
       << bin.lo_ir.transform(0.0) << ", " << bin.lo_ir.transform(1.0)
       << "]\n";

    // Print out number counts
    os << "N:  " << bin.N << " +/- " << bin.N_error
       << " prob: " << bin.N_prob << '\n';
    os << "Nw: " << bin.Nw << " +/- " << bin.Nw_error
       << " prob: " << bin.Nw_prob << '\n';
    os << "Nb: " << bin.Nb << " +/- " << bin.Nb_error
       << " prob: " << bin.Nb_prob << '\n';

    // Print out gamma_ts
    for (auto i = 0u; i < bin.gamma_ts.size(); i++)
      for (auto j = 0u; j < bin.gamma_ts[0].size(); j++)
        os << "gamma_t (src bin " << i << ", R = " << bin.radius[j] << "): "
           << bin.gamma_ts[i][j] << " +/- " << bin.gamma_t_errors[i][j]
           << " prob: " << bin.gamma_t_probs[i][j] << '\n';
  }

  return os;
}

namespace {
  /* Helper function for converting output of integration algorithm into
   * a nicely formatted result.
   */
  template<typename INTEGRAND>
  Gamma_T_Integrated_Bin_Result_S
  make_gamma_t_integrated_bins(const INTEGRAND& gt,
                               const cubacpp::integration_results_v &  results,
                               const size_t n_richness,
                               const size_t n_redshift)
  {
    Gamma_T_Integrated_Bin_Result_S return_arr (n_richness, n_redshift, results.neval, results.nregions, results.status);

    for (auto loi = 0u; loi < n_richness; loi++) {
      for (auto zoi = 0u; zoi < n_redshift; zoi++) {
        return_arr[loi * n_redshift + zoi] = Gamma_T_Integrated_Bin_Result (loi, zoi, gt, results);
      }
    }

    return return_arr;
  }
}

/*
 * The core of this module - A class which represents the integrand for both
 * cluster number counts and gamma_T weak lensing shear.
 *
 * Template parameters:
 *
 * MODELS - A structure representing the specific models to be used in the integration,
 *      see `src/default_models.hh` for an example. Allows easily dropping in your own
 *      model, e.g. if you want to test a different mass-observable relationship,
 *      just subclass DefaultModels with a different MOR type.
 *
 * NRADII - The number of radii at which to compute gamma_T. Ideally this could
 *      be specified at runtime; for now it must be statically known.
 *
 * NRICHNESS, NREDSHIFT - The number of richness and redshift bins to be used,
 *      respectively. As with NRADII, we would like to specify at runtime, but
 *      for now these are compile-time constants. The member functions
 *      `integrate_centered` and `integrate_miscentered` will automatically
 *      sort the integrator output into (richness, redshift) bins, see also
 *      `make_gamma_t_integrated_bins`.
 */
template <typename MODELS>
struct  Gamma_T_Integrand {

  double fcen_;

  typename MODELS::MOR mor;
  typename MODELS::LO_LC lo_lc;
  typename MODELS::LC_LT lc_lt;
  typename MODELS::ZO_ZT zo_zt;
  typename MODELS::PZSOURCE pzsource;
  typename MODELS::ROFFSET roffset;
  typename MODELS::T_CEN T_cen;
  typename MODELS::T_MIS T_mis;
  typename MODELS::A_CEN A_cen;
  typename MODELS::A_MIS A_mis;
  typename MODELS::HMB hmb;
  typename MODELS::HMF hmf;
  typename MODELS::DEL_SIG del_sig;
  typename MODELS::DV_DO_DZ dv_do_dz;
  typename MODELS::OMEGA_Z omega_z;

  y3_cluster::sigma_crit_inv sig_crit_inv_;
  std::shared_ptr<y3_cluster::Interp1D const> da_;

  y3_cluster::IntegrationRange lnM_ir_;
  std::vector<y3_cluster::IntegrationRange> lo_ir_;   /* richness bins */
  // TODO: Matteo uses {(0.3*lo_ir_min), ((2 - loi/8)*lo_ir_max)}
  std::vector<y3_cluster::IntegrationRange> lt_ir_;
  // TODO: Should this be the same as lt_ir_?
  std::vector<y3_cluster::IntegrationRange> lc_ir_;
  std::vector<y3_cluster::IntegrationRange> zo_ir_;   /* redshift bins */
  // TODO: Matteo uses {z_ir_[i].min - 3.5 * z_sigma, z_ir_[i].max + 3.5 * z_sigma}
  std::vector<y3_cluster::IntegrationRange> zt_ir_;
  y3_cluster::IntegrationRange zs_ir_;
  y3_cluster::IntegrationRange R_ir_;
  y3_cluster::IntegrationRange A_ir_;
  y3_cluster::IntegrationRange theta_ir_;

  std::vector<double> r;  /* radii array */
  std::size_t npzsource; /* number of weak lensing source bins */

public:
  // A Gamma_T_Integrand object is constructed by passing in the bunch of
  // callable objects (function pointers or callable class instances)  that
  // specify the various terms of the integrand.
  Gamma_T_Integrand(double fcen,
                    typename MODELS::MOR mor,
                    typename MODELS::LO_LC lo_lc,
                    typename MODELS::LC_LT lc_lt,
                    typename MODELS::ZO_ZT zo_zt,
                    typename MODELS::PZSOURCE pzsource,
                    typename MODELS::ROFFSET roffset,
                    typename MODELS::T_CEN T_cen,
                    typename MODELS::T_MIS T_mis,
                    typename MODELS::A_CEN A_cen,
                    typename MODELS::A_MIS A_mis,
                    typename MODELS::HMB hmb,
                    typename MODELS::HMF hmf,
                    typename MODELS::DEL_SIG del_sig,
                    typename MODELS::DV_DO_DZ dv_do_dz,
                    typename MODELS::OMEGA_Z omega_z,
                    y3_cluster::sigma_crit_inv sci,
                    std::shared_ptr<y3_cluster::Interp1D const> da,
                    y3_cluster::IntegrationRange lnM_ir,
                    std::vector<y3_cluster::IntegrationRange> lo_ir,
                    std::vector<y3_cluster::IntegrationRange> zo_ir,
                    y3_cluster::IntegrationRange zs_ir,
                    y3_cluster::IntegrationRange R_ir,
                    y3_cluster::IntegrationRange A_ir,
                    y3_cluster::IntegrationRange theta_ir,
                    std::vector<double> const& rarray)
    : fcen_(fcen)
    , mor(mor)
    , lo_lc(lo_lc)
    , lc_lt(lc_lt)
    , zo_zt(zo_zt)
    , pzsource(pzsource)
    , roffset(roffset)
    , T_cen(T_cen)
    , T_mis(T_mis)
    , A_cen(A_cen)
    , A_mis(A_mis)
    , hmb(hmb)
    , hmf(hmf)
    , del_sig(del_sig)
    , dv_do_dz(dv_do_dz)
    , omega_z(omega_z)
    , sig_crit_inv_(sci)
    , da_(da)
    , lnM_ir_(lnM_ir)
    , lo_ir_(lo_ir)
    , lt_ir_(y3_cluster::transform(lo_ir, std::function([](y3_cluster::IntegrationRange const& ir) {
                return y3_cluster::IntegrationRange{0.3 * ir.transform(0.0),
                                                    (2 * ir.transform(1.0) > 150)
                                                    ? 150 : 2 * ir.transform(1.0)};
                })))
    , lc_ir_(lt_ir_)
    , zo_ir_(zo_ir)
    , zt_ir_(y3_cluster::transform(zo_ir, std::function([](y3_cluster::IntegrationRange const& ir) {
                const double z_sigma = 0.01;
                return y3_cluster::IntegrationRange{ir.transform(0.0) - 3.5 * z_sigma,
                                                    ir.transform(1.0) + 3.5 * z_sigma};
                })))
    , zs_ir_(zs_ir)
    , R_ir_(R_ir)
    , A_ir_(A_ir)
    , theta_ir_(theta_ir)
    , r(rarray)
    , npzsource(pzsource.nsources())
  {}

  // Alternatively, can automatically construct each model component given a datablock.
  Gamma_T_Integrand(cosmosis::DataBlock& sample,
                    std::vector<double> radii,
                    std::vector<std::shared_ptr<y3_cluster::Interp1D const>> pzsources,
                    std::vector<y3_cluster::IntegrationRange> lo_bins,
                    std::vector<y3_cluster::IntegrationRange> zo_bins)
    : fcen_(get_datablock<double>(sample, "cluster_abundance", "fcen"))
    , mor(sample)
    , lo_lc(sample)
    , lc_lt(sample)
    , zo_zt(sample)
    , pzsource(sample, pzsources)
    , roffset(sample)
    , T_cen(sample)
    , T_mis(sample)
    , A_cen(sample)
    , A_mis(sample)
    , hmb(sample)
    , hmf(sample)
    , del_sig(sample)
    , dv_do_dz(sample)
    , omega_z(sample)
    , sig_crit_inv_(sample)
    , da_(std::make_shared<y3_cluster::Interp1D const>(
                get_datablock<std::vector<double>>(sample, "distances", "z"),
                get_datablock<std::vector<double>>(sample, "distances", "d_a")))
    , lnM_ir_(sample, "lnM")
    , lo_ir_(lo_bins)
    , lt_ir_(y3_cluster::transform(lo_bins, std::function([](y3_cluster::IntegrationRange const& ir) {
                return y3_cluster::IntegrationRange{0.3 * ir.transform(0.0),
                                                    (2 * ir.transform(1.0) > 150)
                                                    ? 150 : 2 * ir.transform(1.0)};
                })))
    , lc_ir_(lt_ir_)
    , zo_ir_(zo_bins)
    , zt_ir_(y3_cluster::transform(zo_bins, std::function([&sample](y3_cluster::IntegrationRange const& ir) {
                const double z_sigma = get_datablock<double>(sample, "cluster_abundance", "zo_zt_sigma");
                return y3_cluster::IntegrationRange{ir.transform(0.0) - 3.5 * z_sigma,
                                                    ir.transform(1.0) + 3.5 * z_sigma};
                })))
    , zs_ir_(sample, "zs")
    , R_ir_(sample, "R")
    , A_ir_(sample, "A")
    , theta_ir_(sample, "theta")
    , r(radii)
    , npzsource(pzsource.nsources())
  {}

  // Convert from one set of bins to another - useful for the
  // `simultaneous_bin_comparison` test and executable, but may prove
  // useful elsewhere
  Gamma_T_Integrand<MODELS>
  with_bins(std::vector<y3_cluster::IntegrationRange> new_lir,
            std::vector<y3_cluster::IntegrationRange> new_zir)
  {
    return {fcen_,
            mor,
            lo_lc,
            lc_lt,
            zo_zt,
            pzsource,
            roffset,
            T_cen,
            T_mis,
            A_cen,
            A_mis,
            hmb,
            hmf,
            del_sig,
            dv_do_dz,
            omega_z,
            sig_crit_inv_,
            da_,
            lnM_ir_,
            // Different lir
            new_lir,
            // Different zir
            new_zir,
            zs_ir_,
            R_ir_,
            A_ir_,
            theta_ir_,
            r};
  }

  typedef std::vector<double> IntegrandResult;

  /* Common integrand functionality. Do not call this directly, you can probably
   * ignore it.
   *
   * Called by both `integrand_centered` and `integrand_miscentered`; factors
   * out calculations common to them both.
   */
  template<typename Fjn, typename Fjg, typename Fnm, typename Fgr>
  IntegrandResult
  integrand_common(double scaled_lt,
                   double scaled_zt,
                   double scaled_zs,
                   double lnM,
                   // Jacobian for N term
                   Fjn jacob_N,
                   // Jacobian for Gamma term
                   Fjg jacob_G,
                   Fnm N_multiplier,
                   // Radially dependent function
                   Fgr gamma_radial_dep) const
  {
    std::size_t const nrichness = lo_ir_.size(),
                      nredshift = zo_ir_.size(),
                      nradii = r.size();

    auto return_arr = IntegrandResult((nradii * npzsource + 3)
                                      * nrichness
                                      * nredshift);

    // We want to minimize the amount of computation - so, avoid recomputation
    // To do this, we will iterate over `z` first and compute only the `z`
    // dependent terms
    for (std::size_t zoi = 0; zoi < nredshift; zoi++) {
      double const zomin = zo_ir_[zoi].transform(0.0);
      double const zomax = zo_ir_[zoi].transform(1.0);
      double const zt = zt_ir_[zoi].transform(scaled_zt);

      const y3_cluster::IntegrationRange this_zs_ir{std::max(zt, zs_ir_.transform(0.0)),
                                                    zs_ir_.transform(1.0)};
      double const zs = this_zs_ir.transform(scaled_zs);

      auto const dv_do_dz_v = dv_do_dz(zt);
      auto const omega_z_v = omega_z(zt);
      auto const hmb_v = hmb(lnM, zt);
      auto const hmf_v = hmf(lnM, zt);
      auto const zo_zt_v = zo_zt(zomin, zomax, zt);

      // These will eventually be passed by CosmoSIS
      double m_shear = 0.0;
      double sig_crit_inv = sig_crit_inv_(zt, zs);
      // This is the lambda-redshift bin weight that we don't fully understand
      double w = 1.0;

      // eq. (28)
      // Compute the z-dependent parts of the gamma_t up front
      //auto gamma_t = std::vector<double>(nradii * npzsource);
      //for (auto i = 0u; i < nradii; i++)
      //  gamma_t[i] = (1.0 + m_shear) / sig_crit * gamma_radial_dep(r[i], zt);
      auto gamma_t = std::vector<double> (nradii * npzsource);
      for (auto srci = 0u; srci < npzsource; srci++) {
        auto const pzsource_v = pzsource(srci, zs);
        for (auto i = 0u; i < nradii; i++) {
          // Nw intentionally left out - returned in return_arr to be used further on
          gamma_t[srci * nradii + i] = (1.0 + m_shear) * pzsource_v * sig_crit_inv * gamma_radial_dep(r[i], zt);
        }
      }

      for (std::size_t loi = 0; loi < nrichness; loi++) {
        auto const richness_bin_start = loi * nredshift * (nradii*npzsource + 3);
        // Zo does not actually need to be integrated over
        double const lt = lt_ir_[loi].transform(scaled_lt);

        auto const mor_v = mor(lt, lnM, zt);

        // eq. (25)
        double const N_int = omega_z_v * dv_do_dz_v * zo_zt_v * hmf_v * mor_v;
        double const N_mult = N_multiplier(loi, zoi);

        // eq. (24)
        double const N = jacob_N(loi, zoi) * N_int * N_mult;
        double const Nw = N * w;
        double const Nb = N * hmb_v;

        // eq. (29)
        auto const gamma_t_int = jacob_G(loi, zoi) * this_zs_ir.jacobian() * N_int * w;

        // eq. (28)
        auto redshift_bin_start = richness_bin_start + zoi * (nradii*npzsource + 3);

        for (auto srci = 0u; srci < npzsource; srci++) {
          for (auto i = 0u; i < nradii; i++) {
            // Nw intentionally left out - returned in return_arr to be used further on
            return_arr[redshift_bin_start + nradii*srci + i] = gamma_t[srci * nradii + i] * gamma_t_int * N_mult;
          }
        }
        return_arr[redshift_bin_start + nradii*npzsource] = N;
        return_arr[redshift_bin_start + nradii*npzsource + 1] = Nw;
        return_arr[redshift_bin_start + nradii*npzsource + 2] = Nb;
      }
    }

    return return_arr;
  }

  /* Miscentered part of integration
   * Term dictionary -
   * * lo - \lambda^{obs}
   * * lc - \lambda^{cen}
   * * lt - \lambda^{true}
   * * zo - z^{obs}
   * * zt - z^{true}
   * * R - R
   * * lnM - ln(M)
   * * A - ???
   * */
  IntegrandResult
  miscentered(double scaled_lo,
              double scaled_lc,
              double scaled_lt,
              double scaled_zt,
              double scaled_zs,
              double scaled_R,
              double scaled_lnM,
              double scaled_A,
              double scaled_theta) const
  {
    // We probably should factor out the common subexpressions, rather than
    // relying upon the optimizer to do a perfect job of this for us. This
    // seems to be the intent of the commented-out code below.
    auto const lnM   = lnM_ir_   .transform(scaled_lnM);
    auto const R     = R_ir_     .transform(scaled_R);
    auto const A     = A_ir_     .transform(scaled_A);
    auto const theta = theta_ir_ .transform(scaled_theta);

    auto jacob_N = [=](std::size_t loi, std::size_t zoi) {
       return lnM_ir_.jacobian() * lo_ir_[loi].jacobian()
              * lt_ir_[loi].jacobian() * lc_ir_[loi].jacobian()
              * zt_ir_[zoi].jacobian()
              * R_ir_.jacobian();
    };
    auto jacob_G = [=](std::size_t loi, std::size_t zoi) {
       return lnM_ir_.jacobian() * lo_ir_[loi].jacobian()
              * lt_ir_[loi].jacobian() * lc_ir_[loi].jacobian()
              * zt_ir_[zoi].jacobian()
              * R_ir_.jacobian() * A_ir_.jacobian()
              * theta_ir_.jacobian();
    };

    // eq. (27)
    auto N_mis = [=](std::size_t loi, std::size_t zoi) {
      auto const lo = lo_ir_[loi].transform(scaled_lo),
                 lc = lc_ir_[loi].transform(scaled_lc),
                 lt = lt_ir_[loi].transform(scaled_lt),
                 zt = zt_ir_[zoi].transform(scaled_zt);
      return (1.0 - fcen_) * lo_lc(lo, lc, R) * lc_lt(lc, lt, zt) * roffset(R);
    };

    // eq. (30)
    // For the following lambda function, `radius` corresponds to what is called
    // `R` in the paper, and `R` corresponds to what is called `R_{mis}` in the
    // paper
    // eq. (31)
    auto gamma_t_mis = [this, N_mis, A, lnM, R, theta](double angle, double zt) {
      // NB: Angle is in arcminutes, must convert to radians
      // TODO: Check this math
      // YZ comment: TODO check unit, this might have been a mpc vs mpc/h thing going on
      double const radius = da_->eval(zt) * angle * pi() / 180.0 / 60.0,
                   adjusted_R = std::sqrt(radius*radius + R*R + 2*R*radius * std::cos(theta));
      return (1.0 / 6.28318530718)
             * exp(A * T_cen(adjusted_R, lnM))/A_ir_.jacobian()
             * del_sig(adjusted_R, lnM, zt);
    };

    return integrand_common(scaled_lt,
                            scaled_zt,
                            scaled_zs,
                            lnM,
                            jacob_N,
                            jacob_G,
                            N_mis,
                            gamma_t_mis);
  }

  /* Centered part of integration
   * Term dictionary -
   * * lo - \lambda^{obs}
   * * lc - \lambda^{cen}
   * * lt - \lambda^{true}
   * * zo - z^{obs}
   * * zt - z^{true}
   * * R - R
   * * lnM - ln(M)
   * * A - ???
   * */
  IntegrandResult
  centered(double scaled_lo,
           double scaled_lt,
           double scaled_zt,
           double scaled_zs,
           double scaled_lnM,
           double scaled_A) const
  {
    // Necessary terms
    auto const lnM = lnM_ir_.transform(scaled_lnM);
    auto const A = A_ir_.transform(scaled_A);

    auto jacob_N = [=](std::size_t loi, std::size_t zoi) {
      return lnM_ir_.jacobian() * lo_ir_[loi].jacobian()
             * lt_ir_[loi].jacobian()
             * zt_ir_[zoi].jacobian();
    };

    auto jacob_G = [=](std::size_t loi, std::size_t zoi) {
      return lnM_ir_.jacobian() * lo_ir_[loi].jacobian()
             * lt_ir_[loi].jacobian()
             * zt_ir_[zoi].jacobian()
             * A_ir_.jacobian();
    };

    // eq. (26)
    auto N_cen = [=](std::size_t loi, std::size_t zoi) {
      auto const lo = lo_ir_[loi].transform(scaled_lo),
                 lt = lt_ir_[loi].transform(scaled_lt),
                 zt = zt_ir_[zoi].transform(scaled_zt);
      return lc_lt(lo, lt, zt) * fcen_;
    };

    // eq. (30)
    // For the following lambda function, `radius` corresponds to what is called
    // `R` in the paper
    auto gamma_t_cen = [this, N_cen, A, lnM](double angle, double zt) {
      // NB: Angle is in arcminutes, must convert to radians
      // TODO: Check this math
      double const radius = da_->eval(zt) * angle * pi() / 180.0 / 60.0;
      return exp(A * T_cen(radius, lnM)) / A_ir_.jacobian() * del_sig(radius, lnM, zt);
    };

    return integrand_common(scaled_lt,
                            scaled_zt,
                            scaled_zs,
                            lnM,
                            jacob_N,
                            jacob_G,
                            N_cen,
                            gamma_t_cen);
  }

  /* Integrates the _centered_ component of the gamma_T, N expressions, and returns
   * a pair of (results, bins), where `results` is the raw `cubacpp` output, and
   * `bins` is an array of `Gamma_T_Integrated_Bin_Result`
   *
   * Arguments:
   *
   * Integrator i: A `cubacpp` integrator, i.e., Cuhre, Vegas, Suave
   * double epsrel: The relative acceptable error of the integration
   * double epsabs: The absolute acceptable integration error
   */
  template<typename Integrator>
  Gamma_T_Integrated_Bin_Result_S
  integrate_centered(Integrator i, double epsrel, double epsabs)
  {
    auto result
        = i.integrate ([this] (double scaled_lo, double scaled_lt,
                               double scaled_zt, double scaled_zs,
                               double scaled_lnM, double scaled_A)
                           {  return centered(scaled_lo,  scaled_lt,
                                              scaled_zt, scaled_zs,
                                              scaled_lnM, scaled_A);  },
                       epsrel, epsabs);

    return make_gamma_t_integrated_bins(*this, result, lo_ir_.size (), zo_ir_.size ());
  }

  /* Integrates the _mis_centered component of the gamma_T, N expressions, and returns
   * a pair of (results, bins), where `results` is the raw `cubacpp` output, and
   * `bins` is an array of `Gamma_T_Integrated_Bin_Result
   *
   * Arguments:
   *
   * Integrator i: A `cubacpp` integrator, i.e., Cuhre, Vegas, Suave
   * double epsrel: The relative acceptable error of the integration
   * double epsabs: The absolute acceptable integration error
   */
  template<typename Integrator>
  Gamma_T_Integrated_Bin_Result_S
  integrate_miscentered(Integrator i, double epsrel, double epsabs)
  {
    auto result = i.integrate([this](double scaled_lo, double scaled_lc, double scaled_lt,
                                     double scaled_zt, double scaled_zs,
                                     double scaled_R, double scaled_lnM,
                                     double scaled_A, double scaled_theta) {
                                 return miscentered(scaled_lo, scaled_lc,
                                                    scaled_lt, scaled_zs,
                                                    scaled_zt, scaled_R, scaled_lnM,
                                                    scaled_A, scaled_theta);
                              },
                              epsrel, epsabs);
    return make_gamma_t_integrated_bins(*this, result, lo_ir_.size (), zo_ir_.size ());
  }
};


template <typename MODELS>
Gamma_T_Integrand<MODELS>
make_gamma_t_integrand(double fcen,
                       typename MODELS::MOR mor,
                       typename MODELS::LO_LC lo_lc,
                       typename MODELS::LC_LT lc_lt,
                       typename MODELS::ZO_ZT zo_zt,
                       typename MODELS::PZSOURCE pzsource,
                       typename MODELS::ROFFSET roffset,
                       typename MODELS::T_CEN t_cen,
                       typename MODELS::T_MIS t_mis,
                       typename MODELS::A_CEN a_cen,
                       typename MODELS::A_MIS a_mis,
                       typename MODELS::HMB hmb,
                       typename MODELS::HMF hmf,
                       typename MODELS::DEL_SIG del_sig,
                       typename MODELS::DV_DO_DZ dv_do_dz,
                       typename MODELS::OMEGA_Z omega_z,
                       std::vector<y3_cluster::IntegrationRange> lo_ir,
                       std::vector<y3_cluster::IntegrationRange> zo_ir,
                       std::size_t  n_radii)
{
  y3_cluster::IntegrationRange lnM_ir{29.0, 38.0};
  y3_cluster::IntegrationRange zt_ir{0.05, 0.35};
  y3_cluster::IntegrationRange zs_ir{0.15, 0.45};
  y3_cluster::IntegrationRange R_ir{0., 3.0};
  y3_cluster::IntegrationRange A_ir{-0.01, 0.01};
  y3_cluster::IntegrationRange theta_ir{0.,6.28318530718};

  auto rarray = std::vector<double> (n_radii);
  for (std::size_t i = 0; i < n_radii; i++)
    rarray[i] = 0.1 * (i + 0.1);

  auto const z_da = read_vector("z.txt"),
             da_arr = read_vector("d_a.txt");

  auto da_f = std::make_shared<y3_cluster::Interp1D const>(z_da, da_arr);

  y3_cluster::sigma_crit_inv sci(da_f);

  return {fcen,
          mor,
          lo_lc,
          lc_lt,
          zo_zt,
          pzsource,
          roffset,
          t_cen,
          t_mis,
          a_cen,
          a_mis,
          hmb,
          hmf,
          del_sig,
          dv_do_dz,
          omega_z,
          sci,
          da_f,
          lnM_ir,
          lo_ir,
          zo_ir,
          zs_ir,
          R_ir,
          A_ir,
          theta_ir,
          rarray };
}


} /*  End of namespace y3_cluster. */


#endif  /* Header guard. */
