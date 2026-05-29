// sigma_prj_gpu_t.cuh -- GPU adaptation of sigma_prj_t.hh (ShearPrjEvaluator)
//
// Costanzi-2026 Sigma_prj / DSigma_prj / gamma_t^prj integrand for CUDA/PAGANI.
// Template-based design matching P_operator_gpu for GPU integrator compatibility.
//
// Integration: 3D over (theta, z, lnM) per (lob_bin, zob, R) grid point.
//
// Template parameter WeightF selects which quantity is integrated:
//   SigmaRnd:  random component of Sigma
//   SigmaCl:   clustering component of Sigma
//   DSigmaRnd: random component of Delta-Sigma
//   DSigmaCl:  clustering component of Delta-Sigma
#ifndef Y3_CLUSTER_CPP_SIGMA_PRJ_GPU_T_CUH
#define Y3_CLUSTER_CPP_SIGMA_PRJ_GPU_T_CUH

#include "utils/make_grid_points.hh"
#include "utils/datablock_reader.hh"

#include "cosmosis/datablock/datablock.hh"
#include "common/cuda/Volume.cuh"

// GPU-compatible models
#include "models/dv_do_dz_t.cuh"
#include "models/hmf_t.cuh"
#include "models/sigma_photoz_des.cuh"

// GPU-compatible interpolators
#include "utils/make_interp_1d.cuh"
#include "utils/make_interp_2d.cuh"

#include <cmath>
#include <optional>
#include <vector>

namespace y3_cuda {

  namespace sp_gpu_detail {

    static constexpr double PI = 3.14159265358979323846;

    __host__ __device__ inline double
    R_lambda(double lam)
    {
      return pow(lam / 100.0, 0.2);
    }

    // DES Y3 richness bin centres: [20,30], [30,45], [45,60], [60,200]
    __host__ __device__ inline double
    lob_center(int bin)
    {
      double const centres[4] = {25.0, 37.5, 52.5, 130.0};
      return (bin >= 0 && bin < 4) ? centres[bin] : 25.0;
    }

    // Photo-z kernel: w_z(z, zob) = 1 - u^2 if |u| < 1, else 0
    __host__ __device__ inline double
    w_photoz(double z, double zob, double sigma_z)
    {
      double const u = (z - zob) / sigma_z;
      if (fabs(u) >= 1.0) return 0.0;
      return 1.0 - u * u;
    }

    // Per-z LoS-slab exclusion angle
    __host__ __device__ inline double
    theta_excl_at_z(double chi_z, double chi_o, double R_excl)
    {
      double const denom = 2.0 * chi_z * chi_o + 1.0e-30;
      double cos_ex = (chi_z * chi_z + chi_o * chi_o - R_excl * R_excl) / denom;
      if (cos_ex >= 1.0 - 1.0e-12) return 0.0;
      if (cos_ex < -1.0) cos_ex = -1.0;
      return acos(cos_ex);
    }

    // Sigmoid for b_sel(theta) = B_small + (B_large - B_small) * sigmoid
    __host__ __device__ inline double
    b_sel_sigmoid(double theta, double theta_lob)
    {
      double const k = 2.5 / theta_lob;
      double const theta0 = 0.5 * theta_lob;
      return 1.0 / (1.0 + exp(-k * (theta - theta0)));
    }

    // Simplified NFW Sigma (placeholder - real version uses lookup tables)
    __host__ __device__ inline double
    compute_sigma_nfw(double R, double R_mis, double lnM)
    {
      double const M = exp(lnM);
      double const r_s = 0.1 * pow(M / 1e14, 0.33);
      double const x = sqrt(R * R + R_mis * R_mis) / r_s;
      if (x < 0.01) return 1.0e15 / (x + 0.01);
      if (x > 100.0) return 0.0;
      return 1.0e14 * pow(M / 1e14, 0.8) / (x * (1.0 + x) * (1.0 + x));
    }

    __host__ __device__ inline double
    compute_dsigma_nfw(double R, double R_mis, double lnM)
    {
      double const M = exp(lnM);
      double const r_s = 0.1 * pow(M / 1e14, 0.33);
      double const x = sqrt(R * R + R_mis * R_mis) / r_s;
      if (x < 0.01) return 1.0e15 / (x + 0.01);
      if (x > 100.0) return 0.0;
      return 1.0e14 * pow(M / 1e14, 0.8) / (x * x * (1.0 + x));
    }

  }  // namespace sp_gpu_detail


  // =====================================================================
  // GPU integrand template for Sigma_prj / DSigma_prj
  // =====================================================================
  template <typename WeightF>
  class SigmaPrj_gpu {
  public:
    using grid_t       = y3_cluster::grid_t<4>;  // (lambda_bin, zo_low, zo_high, radii)
    using grid_point_t = grid_t::value_type;

  private:
    using volume_t = quad::Volume<double, 3>;    // (theta, z, lnM)

    // GPU-compatible models
    std::optional<y3_cuda::HMF_t>        hmf_;
    std::optional<y3_cuda::DV_DO_DZ_t>   dv_do_dz_;

    // GPU-compatible interpolators
    std::optional<quad::Interp2D> hmb_;       // haloModel/bias(lnM, z)
    std::optional<quad::Interp2D> xi_nl_;     // xi_nl(r, z)
    std::optional<quad::Interp1D> chi_;       // distances/d_c(z)
    std::optional<quad::Interp1D> sci_;       // average_sigma_crit_inv(z)

    // b_sel asymptotes from bsel.py
    std::optional<quad::Interp2D> b_small_;   // b_sel_marginalised/b_small(lob, zob)
    std::optional<quad::Interp2D> b_large_;   // b_sel_marginalised/b_large(lob, zob)

    // CPU-side sci interpolation data (for get_sci() called from host)
    std::vector<double> sci_z_;
    std::vector<double> sci_vals_;

    // Current grid point state
    double zob_         = 0.0;
    double R_           = 0.0;
    int    lob_bin_     = 0;

    // Cosmology
    double h0_ = 1.0;

  public:
    explicit SigmaPrj_gpu(cosmosis::DataBlock& /*cfg*/)
    {
      // No special setup needed
    }

    void
    set_sample(cosmosis::DataBlock& sample)
    {
      hmf_.emplace(sample);
      dv_do_dz_.emplace(sample);

      // Interpolators from datablock
      hmb_.emplace(make_Interp2D(sample, "haloModel", "lnM", "z", "bias"));
      xi_nl_.emplace(make_Interp2D(sample, "xi_nl", "r", "z", "xi_nl"));
      chi_.emplace(make_Interp1D(sample, "distances", "z", "d_c"));

      // Optional average_sigma_crit_inv for gamma_t (CPU-side for get_sci())
      if (sample.has_val("average_sigma_crit_inv", "zlense")) {
        sci_.emplace(make_Interp1D(sample, "average_sigma_crit_inv",
                                   "zlense", "sci_average"));
        // Also store CPU-side copy for get_sci() called from host
        sci_z_ = sample.view<std::vector<double>>("average_sigma_crit_inv", "zlense");
        sci_vals_ = sample.view<std::vector<double>>("average_sigma_crit_inv", "sci_average");
      }

      // b_sel asymptotes from bsel.py (2D: shape is (n_zob, n_lob))
      if (sample.has_val("b_sel_marginalised", "b_small") &&
          sample.has_val("b_sel_marginalised", "b_large")) {
        auto lob_vec = sample.view<std::vector<double>>("b_sel_marginalised", "lob");
        auto zob_vec = sample.view<std::vector<double>>("b_sel_marginalised", "zob");
        auto const& b_small_nd = sample.view<cosmosis::ndarray<double>>(
            "b_sel_marginalised", "b_small");
        auto const& b_large_nd = sample.view<cosmosis::ndarray<double>>(
            "b_sel_marginalised", "b_large");

        std::vector<double> b_small_vec(b_small_nd.begin(), b_small_nd.end());
        std::vector<double> b_large_vec(b_large_nd.begin(), b_large_nd.end());

        b_small_.emplace(quad::Interp2D(zob_vec, lob_vec, b_small_vec));
        b_large_.emplace(quad::Interp2D(zob_vec, lob_vec, b_large_vec));
      }

      h0_ = sample.view<double>("cosmological_parameters", "h0");
    }

    void
    set_grid_point(grid_point_t const& pt)
    {
      // Only store values - do NOT call interpolators here (they run on GPU)
      lob_bin_  = static_cast<int>(pt[0]);
      zob_      = 0.5 * (pt[1] + pt[2]);
      R_        = pt[3];
    }

    size_t
    get_device_mem_footprint()
    {
      size_t size = 0;
      if (hmf_) size += hmf_->get_device_mem_footprint();
      if (dv_do_dz_) size += dv_do_dz_->get_device_mem_footprint();
      return size;
    }

    // Integrand at (theta, z, lnM) -- returns scalar based on WeightF
    __host__ __device__ double
    operator()(double theta, double z, double lnM) const
    {
      using namespace sp_gpu_detail;

      // Photo-z kernel
      SIGMA_PHOTOZ_DES_t sigma_z_model;
      double const sig_z = sigma_z_model(z);
      double const wz = w_photoz(z, zob_, sig_z);
      if (wz <= 0.0) return 0.0;

      // Compute geometry (moved from set_grid_point - these use GPU interpolators)
      double const lobc = lob_center(lob_bin_);
      double const chi_o = chi_->clamp(zob_) * h0_;
      double const D_A_o = chi_o / (1.0 + zob_);
      double const R_excl = R_lambda(lobc) * (1.0 + zob_);
      double const theta_lob = R_lambda(lobc) * (1.0 + zob_) / chi_o;

      // b_sel asymptotes (use fallback if not available)
      double B_small = 1.0;
      double B_large = 1.0;
      if (b_small_ && b_large_) {
        B_small = b_small_->clamp(zob_, lobc);
        B_large = b_large_->clamp(zob_, lobc);
      }

      // Geometry for this z
      double const chi_z = chi_->clamp(z) * h0_;
      double const theta_excl_z = theta_excl_at_z(chi_z, chi_o, R_excl);
      double const sin_th = sin(theta);
      double const cos_th = cos(theta);

      // Common weights
      double const dV = (*dv_do_dz_)(z);
      double const hmf_val = (*hmf_)(lnM, z);
      double const common = 2.0 * PI * sin_th * dV * wz * hmf_val;
      if (common <= 0.0) return 0.0;

      // NFW Sigma_mis / DSigma_mis
      double const R_th = theta * D_A_o;
      double const Sigma_v = compute_sigma_nfw(R_, R_th, lnM);
      double const DSigma_v = compute_dsigma_nfw(R_, R_th, lnM);

      // Clustering component (only if theta > theta_excl)
      double xi_val = 0.0;
      double b_Mz = 0.0;
      double b_sel = 0.0;
      if (theta > theta_excl_z) {
        double const dchi = sqrt(fmax(
            chi_z * chi_z + chi_o * chi_o
            - 2.0 * chi_z * chi_o * cos_th, 0.0));
        xi_val = xi_nl_->clamp(dchi, zob_);
        b_Mz = hmb_->clamp(lnM, z);
        b_sel = B_small + (B_large - B_small)
                         * b_sel_sigmoid(theta, theta_lob);
      }

      // Return the selected component
      return WeightF::weight(common, Sigma_v, DSigma_v, b_Mz, b_sel, xi_val,
                             theta > theta_excl_z);
    }

    // CPU-side linear interpolation for sci (called from host code)
    double get_sci() const {
      if (sci_z_.empty()) return 0.0;
      // Simple linear interpolation
      if (zob_ <= sci_z_.front()) return sci_vals_.front();
      if (zob_ >= sci_z_.back()) return sci_vals_.back();
      for (size_t i = 1; i < sci_z_.size(); ++i) {
        if (zob_ <= sci_z_[i]) {
          double const t = (zob_ - sci_z_[i-1]) / (sci_z_[i] - sci_z_[i-1]);
          return sci_vals_[i-1] + t * (sci_vals_[i] - sci_vals_[i-1]);
        }
      }
      return sci_vals_.back();
    }

    static char const* module_label() { return "ShearPrjGPU"; }

    // Create 3D integration volumes from ini config
    static std::vector<volume_t>
    make_integration_volumes(cosmosis::DataBlock& cfg)
    {
      char const* label = module_label();

      double theta_low  = cfg.view<double>(label, "theta_low");
      double theta_high = cfg.view<double>(label, "theta_high");
      double zt_low     = cfg.view<double>(label, "zt_low");
      double zt_high    = cfg.view<double>(label, "zt_high");
      double lnm_low    = cfg.view<double>(label, "lnm_low");
      double lnm_high   = cfg.view<double>(label, "lnm_high");

      volume_t vol;
      vol.lows[0]  = theta_low;   vol.highs[0] = theta_high;
      vol.lows[1]  = zt_low;      vol.highs[1] = zt_high;
      vol.lows[2]  = lnm_low;     vol.highs[2] = lnm_high;

      return {vol};
    }

    // Create grid points from ini config
    static grid_t
    make_grid_points(cosmosis::DataBlock& cfg)
    {
      char const* label = module_label();

      auto lambda_bin = get_vector_double(cfg, label, "lambda_bin");
      auto zo_low     = get_vector_double(cfg, label, "zo_low");
      auto zo_high    = get_vector_double(cfg, label, "zo_high");
      auto radii      = get_vector_double(cfg, label, "radii");

      size_t n = lambda_bin.size();
      if (zo_low.size() != n || zo_high.size() != n || radii.size() != n) {
        throw std::runtime_error(
          std::string(label) + ": grid arrays must have same length");
      }

      grid_t result;
      result.set_names({"lambda_bin", "zo_low", "zo_high", "radii"});
      for (size_t i = 0; i < n; ++i) {
        result.push_back({lambda_bin[i], zo_low[i], zo_high[i], radii[i]});
      }
      return result;
    }
  };

  // ---------------------------------------------------------------------
  // Weight functors for GPU
  // ---------------------------------------------------------------------

  struct SigmaRndWeight {
    static char const* name() { return "sigma_rnd"; }
    __host__ __device__ static double
    weight(double common, double Sigma_v, double /*DSigma_v*/,
           double /*b_Mz*/, double /*b_sel*/, double /*xi_val*/, bool /*has_cl*/)
    {
      return common * Sigma_v;
    }
  };

  struct SigmaClWeight {
    static char const* name() { return "sigma_cl"; }
    __host__ __device__ static double
    weight(double common, double Sigma_v, double /*DSigma_v*/,
           double b_Mz, double b_sel, double xi_val, bool has_cl)
    {
      if (!has_cl) return 0.0;
      return common * b_Mz * b_sel * xi_val * Sigma_v;
    }
  };

  struct DSigmaRndWeight {
    static char const* name() { return "dsigma_rnd"; }
    __host__ __device__ static double
    weight(double common, double /*Sigma_v*/, double DSigma_v,
           double /*b_Mz*/, double /*b_sel*/, double /*xi_val*/, bool /*has_cl*/)
    {
      return common * DSigma_v;
    }
  };

  struct DSigmaClWeight {
    static char const* name() { return "dsigma_cl"; }
    __host__ __device__ static double
    weight(double common, double /*Sigma_v*/, double DSigma_v,
           double b_Mz, double b_sel, double xi_val, bool has_cl)
    {
      if (!has_cl) return 0.0;
      return common * b_Mz * b_sel * xi_val * DSigma_v;
    }
  };

}  // namespace y3_cuda

#endif
