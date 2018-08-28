#ifndef Y3_CLUSTER_SIGMA_HH
#define Y3_CLUSTER_SIGMA_HH

#include <cmath>
#include <cubacpp/cubacpp.hh>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <memory>
#include <stdexcept>
#include <vector>

#include "bessel_polynomial_integral.hh"
#include "datablock_reader.hh"
#include "integration_range.hh"
#include "interp_1d.hh"
#include "interp_2d.hh"
#include "omega_z_sdss.hh"
#include "read_vector.hh"
#include "primitives.hh"

namespace y3_cluster {
  // Coefficients of survey mask, in spherical harmonics
  double
  survey_mask_kay(unsigned l, double z)
  {
    OMEGA_Z_SDSS omega_z;

    if (l == 0)
        return 0.5 / std::sqrt(pi());
    else {
      const double oz = omega_z(z),
                   cos_theta = 1 - oz / 2 / pi(),
                   frac = (std::legendre(l - 1, cos_theta) - std::legendre(l + 1, cos_theta))
                          / oz;
      return std::sqrt(pi() / (2*l + 1)) * frac;
    }
  }

  namespace {
    // Create a grid of K_l(z_i), so: compute_ks()[i][j] = K_j(z_i)
    std::vector<std::vector<double>>
    compute_ks(unsigned maxl, const std::vector<IntegrationRange>& zirs)
    {
      std::vector<std::vector<double>> out;
      for (const auto zbin : zirs) {
        std::vector<double> this_bin;
        // Use the midpoint of the bin
        for (auto l = 0u; l < maxl; l++)
          this_bin.push_back(survey_mask_kay(l, zbin.transform(0.5)));
        out.push_back(this_bin);
      }
      return out;
    }
  }

  /* The Sample Variance covariance is computed by:
   *
   *  C_{ij} = Nb_i * Nb_j * sigma_{ij}^2
   *
   * The class in this file computes a matrix of the last term - sigma_{ij}^2
   */
  class SampleVariance_t {
    std::vector<IntegrationRange> z_ranges;
    std::shared_ptr<Interp2D const> matter_power_lin;
    std::shared_ptr<Interp1D const> dcom;

    // Number of Bessel functions to sum over, i.e. l = 0, 1, ..., maxl
    std::size_t maxl;
    std::vector<std::vector<double>> ks;

    // Reduced hubble parameter `h`
    double hubble;

  public:
    SampleVariance_t() = delete;
    SampleVariance_t(const std::vector<IntegrationRange>& ir, double h)
      : z_ranges(ir)
      , matter_power_lin(std::make_shared<Interp2D const>(
            read_vector("matter_power_lin/k_h.txt"),
            read_vector("matter_power_lin/z.txt"),
            read_vector("matter_power_lin/p_k.txt")))
      , dcom(std::make_shared<Interp1D const>(
            read_vector("distances/z.txt"),
            read_vector("distances/d_m.txt")))
      , maxl(5)
      , ks(compute_ks(maxl, ir))
      , hubble(h)
      {}

    // TODO: Implement grabbing z-bins from datablock, and maxl parameter
    explicit SampleVariance_t(cosmosis::DataBlock& sample)
      : z_ranges()
      , matter_power_lin(std::make_shared<Interp2D const>(
            get_datablock<std::vector<double>>(sample, "matter_power_lin", "k_h"),
            get_datablock<std::vector<double>>(sample, "matter_power_lin", "z"),
            get_datablock<cosmosis::ndarray<double>>(sample, "matter_power_lin", "p_k")))
      , dcom(std::make_shared<Interp1D const>(
            get_datablock<std::vector<double>>(sample, "distances", "z"),
            get_datablock<std::vector<double>>(sample, "distances", "d_m")))
      , maxl(5)
      , ks(compute_ks(maxl, {}))
      , hubble(get_datablock<double>(sample, "cosmological_parameters", "h0"))
      {}

    /* Computes:
     *
     *    sum_{l=0}^{maxl} K_l(i) * K_l(j) * R_i(k) * R_j(k)
     */
    double
    compute_sum_over_bessels(double k, unsigned i, unsigned j) const
    {
      const double ri_min = dcom->eval(z_ranges[i].transform(0)) * hubble,
                   ri_max = dcom->eval(z_ranges[i].transform(1)) * hubble,
                   rj_min = dcom->eval(z_ranges[j].transform(0)) * hubble,
                   rj_max = dcom->eval(z_ranges[j].transform(1)) * hubble,
                   i_normalization = (ri_max * ri_max * ri_max - ri_min * ri_min * ri_min) / 3.0,
                   j_normalization = (rj_max * rj_max * rj_max - rj_min * rj_min * rj_min) / 3.0;

      const std::vector<double> bessels_i = bessel_polynomial_integrals(2, maxl, k, ri_min, ri_max),
                                bessels_j = (i == j) ? bessels_i :
                                    bessel_polynomial_integrals(2, maxl, k, ri_min, ri_max);

      double sum = 0;
      for (auto l = 0u; l < maxl; l++) {
        double this_l = bessels_i[l] * bessels_j[l]
                      / i_normalization / j_normalization
                      * ks[i][l] * ks[j][l];
        if (l & 1)
          sum -= this_l;
        else
          sum += this_l;
      }
      return sum;
    }

    double
    manual_compute_sum_over_bessels(double k, unsigned i, unsigned j) const
    {
      const double ri_min = dcom->eval(z_ranges[i].transform(0)) * hubble,
                   ri_max = dcom->eval(z_ranges[i].transform(1)) * hubble,
                   rj_min = dcom->eval(z_ranges[j].transform(0)) * hubble,
                   rj_max = dcom->eval(z_ranges[j].transform(1)) * hubble,
                   i_normalization = (ri_max * ri_max * ri_max - ri_min * ri_min * ri_min) / 3.0,
                   j_normalization = (rj_max * rj_max * rj_max - rj_min * rj_min * rj_min) / 3.0;

      cubacpp::Cuhre integrator;
      integrator.maxeval = 999999999;
      auto integrate_bessel = [&](unsigned l, double k, double r_min, double r_max) {
          y3_cluster::IntegrationRange r_ir{r_min, r_max};
          const auto result = integrator.integrate([&](double r_scaled, double ignored [[maybe_unused]]) {
                  const double r = r_ir.transform(r_scaled);
                  return r * r * gsl_sf_bessel_jl(l, k * r) * r_ir.jacobian();
                  }, 1e-2, 1e-18);

          if (result.status != 0)
              throw std::runtime_error("Big boom! (bessel integral did not converge)");

          return result.value;
      };

      double sum = 0;
      for (auto l = 0u; l < maxl; l++) {
        double this_l = integrate_bessel(l, k, ri_min, ri_max) * integrate_bessel(l, k, rj_min, rj_max)
                      / i_normalization / j_normalization
                      * ks[i][l] * ks[j][l];
        if (l & 1)
          sum -= this_l;
        else
          sum += this_l;
      }
      return sum;
    }

    // Compute the full sigma_{ij} matrix
    std::vector<std::vector<double>>
    compute() const
    {
      std::vector<std::vector<double>> sigma_mat(z_ranges.size());
      for (auto& v : sigma_mat)
        v.resize(z_ranges.size());

      // Compute matrix of sigma^{SampVar}_{ij}
      cubacpp::Cuhre integrator;
      integrator.maxeval = 999999999;

      // Using Matteo's range for now
      IntegrationRange lnk_ir{std::log(0.0001), std::log(0.80)};
      for (auto i = 0u; i < z_ranges.size(); i++) {
        for (auto j = 0u; j < z_ranges.size(); j++) {
          if (i <= j) {
            /* Integrate \sigma_{ij} =
             *      ((dV/d\Omega)|_{\Delta zi} (dV/d\Omega)|_{\Delta zi})^{-1}
             *       2/pi \int dk^3 \int_{\Delta z_i} dz_i \int_{\Delta z_j} dz_j
             *       \sqrt{P_{lin}(k, z_i) P_{lin}(k, z_j)} D(z_i)^2 D(z_j)^2
             *       \sum_l (-1}^l j_l(k D(z_i)) j_l(k D(z_j)) K_l(z_i) K_l(z_j)
             */

            // Have ignored parameter so we can use Cuhre - which requires >=2 integration dimensions
            // But, Cuhre is much, much faster at this than Vegas, Suave
            auto result = integrator.integrate([&](const double lnk_scaled, const double ignored [[maybe_unused]]) {
                      const auto lnk = lnk_ir.transform(lnk_scaled),
                                 zi_mid = z_ranges[i].transform(0.5),
                                 zj_mid = z_ranges[j].transform(0.5),
                                 k = std::exp(lnk);

                      const auto matter_power = std::sqrt(matter_power_lin->eval(k, zi_mid) * matter_power_lin->eval(k, zj_mid)),
                                 sum = compute_sum_over_bessels(k, i, j);

                      // k^3 due to integration over log
                      return 2 / pi()
                             * k * k * k
                             * sum * matter_power
                             * lnk_ir.jacobian();
                    },
                    1e-2, 1e-18);

            if (result.status != 0)
              throw std::runtime_error("Integration of dvdodz failed in SampleVariance");

            sigma_mat[i][j] = result.value;
          } else {
            // Already computed (sigma_{ij} = sigma_{ji})
            sigma_mat[i][j] = sigma_mat[j][i];
          }
        }
      }

      return sigma_mat;
    }
  };
}

#endif
