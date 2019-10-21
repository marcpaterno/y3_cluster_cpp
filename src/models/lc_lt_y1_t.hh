#ifndef Y3_CLUSTER_LC_LT_Y1_T_HH
#define Y3_CLUSTER_LC_LT_Y1_T_HH

#include "cosmosis/datablock/datablock.hh"
#include "utils/interp_2d.hh"
#include "utils/primitives.hh"
#include "utils/read_vector.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

namespace y3_cluster {
  struct LC_LT_Y1_t {

    std::shared_ptr<Interp2D const> tau_interp;
    std::shared_ptr<Interp2D const> mu_interp;
    std::shared_ptr<Interp2D const> sigma_interp;
    std::shared_ptr<Interp2D const> fmsk_interp;
    std::shared_ptr<Interp2D const> fprj_interp;

    LC_LT_Y1_t();
    // LC_LT_Y1_t(const cosmosis::DataBlock&) {TODO}

    double
    operator()(double lamcent, double lamtrue, double ztrue) const
    {
      // Values of each fit parameter at lamtrue, ztrue
      const auto tau = tau_interp->eval(lamtrue, ztrue);
      const auto mu = mu_interp->eval(lamtrue, ztrue);
      const auto sigma = sigma_interp->eval(lamtrue, ztrue);
      const auto fmsk = fmsk_interp->eval(lamtrue, ztrue);
      const auto fprj = fprj_interp->eval(lamtrue, ztrue);

      // Some repeated expressions
      const auto exptau =
        std::exp(tau * (2.0 * mu + tau * sigma * sigma - 2.0 * lamcent) / 2.0);
      const auto root_two_sigma = std::sqrt(2.0) * sigma;
      const auto mu_tau_sig_sqr = mu + tau * sigma * sigma;

      // Helper function for common pattern
      const auto erfc_scaled = [root_two_sigma](double a, double b) {
        return std::erfc((a - b) / root_two_sigma);
      };

      return (1.0 - fmsk) * (1.0 - fprj) * y3_cluster::gaussian(lamcent, mu, sigma) +
             0.5 * ((1.0 - fmsk) * fprj * tau + fmsk * fprj / lamtrue) * exptau *
               erfc_scaled(mu_tau_sig_sqr, lamcent) +
             0.5 * fmsk / lamtrue *
               (erfc_scaled(lamcent, mu) - erfc_scaled(lamcent + lamtrue, mu)) -
             0.5 * fmsk * fprj / lamtrue *
               (std::exp(-tau * lamtrue) * exptau *
                erfc_scaled(mu_tau_sig_sqr, lamcent + lamtrue));
    }
  };


  LC_LT_Y1_t::LC_LT_Y1_t() {
    // This should be read in from the data block
    std::string runname = "DESY1A_v1.4";
    std::string fileroot = "lc_lt_fits/" + runname;
  
    auto lamtrue_in = read_vector(fileroot + "/lamtrue.txt");
    auto z_in = read_vector(fileroot + "/z.txt");

    // auto tau_in = read_vector(fileroot + "tau.txt");
    // auto mu_in = read_vector(fileroot + "mu.txt");
    // auto sigma_in = read_vector(fileroot + "sigma.txt");
    // auto fmsk_in = read_vector(fileroot + "fmsk.txt");
    // auto fprj_in = read_vector(fileroot + "fprj.txt");

    LC_LT_Y1_t::tau_interp = std::make_shared<Interp2D const>(
        lamtrue_in, z_in, read_vector(fileroot + "/tau.txt"));
    LC_LT_Y1_t::mu_interp = std::make_shared<Interp2D const>(
        lamtrue_in, z_in, read_vector(fileroot + "/mu.txt"));
    LC_LT_Y1_t::sigma_interp = std::make_shared<Interp2D const>(
        lamtrue_in, z_in, read_vector(fileroot + "/sigma.txt"));
    LC_LT_Y1_t::fmsk_interp = std::make_shared<Interp2D const>(
        lamtrue_in, z_in, read_vector(fileroot + "/fmsk.txt"));
    LC_LT_Y1_t::fprj_interp = std::make_shared<Interp2D const>(
        lamtrue_in, z_in, read_vector(fileroot + "/fprj.txt"));
  }

}









#endif
