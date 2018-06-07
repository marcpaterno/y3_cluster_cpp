#include "test/lc_lt_t.hh"

using namespace y3_cluster;

namespace {
  std::vector<double> const xxs(LC_LT_t::lt_bins.begin(),
                                LC_LT_t::lt_bins.end() - 1);
  std::vector<double> const yys(LC_LT_t::zt_bins.begin(),
                                LC_LT_t::zt_bins.end() - 1);
  std::vector<double> const taus(LC_LT_t::tau_arr.begin(),
                                 LC_LT_t::tau_arr.end());
  std::vector<double> const mus(LC_LT_t::mu_arr.begin(), LC_LT_t::mu_arr.end());
  std::vector<double> const sigmas(LC_LT_t::sigma_arr.begin(),
                                   LC_LT_t::sigma_arr.end());
  std::vector<double> const fmsks(LC_LT_t::fmsk_arr.begin(),
                                  LC_LT_t::fmsk_arr.end());
  std::vector<double> const fprjs(LC_LT_t::fprj_arr.begin(),
                                  LC_LT_t::fprj_arr.end());

}

Interp2D const LC_LT_t::tau_interp(xxs, yys, taus);
Interp2D const LC_LT_t::mu_interp(xxs, yys, mus);
Interp2D const LC_LT_t::sigma_interp(xxs, yys, sigmas);
Interp2D const LC_LT_t::fmsk_interp(xxs, yys, fmsks);
Interp2D const LC_LT_t::fprj_interp(xxs, yys, fprjs);
