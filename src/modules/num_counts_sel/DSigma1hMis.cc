// N_i[DeltaSigma_cl(R)] -- 1-halo differential surface density with
// miscentering (DES-Y3 redMaPPer):
//
//   DSigma_cl(R, M, z) = (1 - f_mis) * DSigma_NFW(R, M)
//                      + f_mis      * DSigma_mis(R, M; tau_mis * R_lambda)
//
// Ratio DSigma1hMisSel / NumCountsSel gives the bin-averaged
// <DeltaSigma_cl(R)> including the miscentering tail.  Mirrors the
// centred DSigma1hSel driver; the only difference is the weight
// functor, which carries the gamma-kernel NFW_DSIGMA_MIS table and
// the (f_mis, tau_mis) mixture parameters.
#include "models/n_operator_sel_t.hh"
#include "modules/num_counts_sel/lensing_weights.hh"
#include "utils/module_macros.hh"

using DSigma1hMisSel =
  y3_cluster::NOperatorSelRadial<y3_cluster_sel_weights::DSigma1hMisWeight>;

DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE(DSigma1hMisSel)
