#ifndef Y3_CLUSTER_SIGMA_HH
#define Y3_CLUSTER_SIGMA_HH

#include <cmath>
#include <cubacpp/cubacpp.hh>
#include <memory>
#include <stdexcept>
#include <vector>

#include "datablock_reader.hh"
#include "dv_do_dz_t.hh"
#include "ez.hh"
#include "integration_range.hh"
#include "interp_2d.hh"
#include "omega_z_sdss.hh"
#include "primitives.hh"

namespace y3_cluster {
  namespace {
    std::vector<double>
    convert_to_mpch(std::vector<double> dm, double h)
    {
      std::vector<double> out(dm);
      for (auto& v : out)
          v *= h;
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
    DV_DO_DZ_t dvdodz;
    OMEGA_Z_SDSS omega_z;
    EZ ez;
    std::shared_ptr<Interp2D const> matter_power_lin;
    std::shared_ptr<Interp1D const> dcom;

    std::size_t maxl;

  public:
    SampleVariance_t() = delete;
    SampleVariance_t(EZ ez, DV_DO_DZ_t dvdodz, const std::vector<IntegrationRange>& ir)
      : z_ranges(ir)
      , dvdodz(dvdodz)
      , omega_z()
      , ez(ez)
      , matter_power_lin(std::make_shared<Interp2D const>(
            read_vector("matter_power_lin/k_h.txt"),
            read_vector("matter_power_lin/z.txt"),
            read_vector("matter_power_lin/p_k.txt")))
      , dcom(std::make_shared<Interp1D const>(
            read_vector("distances/z.txt"),
            read_vector("distances/d_m.txt")))
      , maxl(100)
      {}

    explicit SampleVariance_t(cosmosis::DataBlock& sample)
      : z_ranges()
      , dvdodz(sample)
      , omega_z()
      , ez(sample)
      , matter_power_lin(std::make_shared<Interp2D const>(
            get_datablock<std::vector<double>>(sample, "matter_power_lin", "k_h"),
            get_datablock<std::vector<double>>(sample, "matter_power_lin", "z"),
            get_datablock<cosmosis::ndarray<double>>(sample, "matter_power_lin", "p_k")))
      , dcom(std::make_shared<Interp1D const>(
            get_datablock<std::vector<double>>(sample, "distances", "z"),
            convert_to_mpch(get_datablock<std::vector<double>>(sample, "distances", "d_m"),
                            get_datablock<double>(sample, "distances", "h"))))
      , maxl(100)
      {}

    double
    kay(std::size_t l, double z) const
    {
      if (l == 0)
          return 1 / 2 / std::sqrt(pi());
      else {
        const double oz = omega_z(z),
                     cos_theta = 1 - oz / 2 / pi(),
                     frac = (std::legendre(l - 1, cos_theta) - std::legendre(l + 1, cos_theta))
                            / oz;
        return std::sqrt(pi() / (2*l + 1)) * frac;
      }
    }

    /* Compute the full sigma_{ij} matrix
     */
    std::vector<std::vector<double>>
    operator()()
    {
      std::vector<std::vector<double>> sigma_mat(z_ranges.size());
      for (auto& v : sigma_mat)
        v.resize(z_ranges.size());

      /* Precompute:
       *
       *  (dV/d\Omega)|_{\Delta z_i} = \int_{\Delta z_i} dz (dV/d\Omega/dz)(z)
       */
      std::vector<double> delta_dvdos(z_ranges.size());
      cubacpp::QAG qag;
      for (auto i = 0u; i < z_ranges.size(); i++) {
        auto result = qag.with_range(z_ranges[i].transform(0.0),
                                     z_ranges[i].transform(1.0))
                         .integrate([this] (double z) { return dvdodz(z); },
                                    1e-5, 1e-18);

        if (result.status != 0)
          throw std::runtime_error("Integration of dvdodz failed in SampleVariance");

        delta_dvdos[i] = result.value;
      }


      /* Compute matrix of sigma^{SampVar}_{ij} */
      cubacpp::Cuhre cuhre;
      cuhre.maxeval = 999999999;

      // No idea what real range should be
      IntegrationRange lnk_ir{1, 10};
      for (auto i = 0u; i < z_ranges.size(); i++) {
        for (auto j = 0u; j < z_ranges.size(); j++) {
          if (i <= j) {
            // Integrate \sigma_{ij} =
            // TODO: Implement dcom, bessel, kay, and matter_power_lin
            auto result = cuhre.integrate([&](double lnk_scaled, double zi_scaled, double zj_scaled) {
                      const auto lnk = lnk_ir.transform(lnk_scaled),
                                 zi = z_ranges[i].transform(zi_scaled),
                                 zj = z_ranges[j].transform(zj_scaled),
                                 k = std::exp(lnk);

                      const auto matter_power = std::sqrt(matter_power_lin->eval(k, zi) * matter_power_lin->eval(k, zj)),
                                 dcomi = dcom->eval(zi),
                                 dcomj = dcom->eval(zj);

                      double sum = 0;
                      for (auto l = 0u; l < maxl; l++) {
                        const double this_l = std::sph_bessel(l, k * dcomi) * std::sph_bessel(l, k * dcomj)
                                            * kay(l, zi) * kay(l, zj);
                        if (l & 1)
                          sum -= this_l;
                        else
                          sum += this_l;
                      }

                      // k^3 due to integration over log
                      return k * k * k * dcomi * dcomi * dcomj * dcomj
                             * sum * matter_power
                             * lnk_ir.jacobian() * z_ranges[i].jacobian() * z_ranges[j].jacobian();
                    },
                    1e-3, 1e-18);

            if (result.status != 0)
              throw std::runtime_error("Integration of dvdodz failed in SampleVariance");

            sigma_mat[i][j] = 2 * result.value / delta_dvdos[i] / delta_dvdos[j] / pi();
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
