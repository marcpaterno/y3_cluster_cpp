// N_i[gamma_t^1h, full(R)] -- 1-halo tangential shear with miscentering:
//
//   gamma_t^1h, full(R, M, z) = DSigma_cl(R, M, z) * Sigma_crit^-1(z_lens)
//
// where DSigma_cl is the (centred + miscentered) mixture defined in
// DSigma1hMis.cc.  Ratio Shear1hMisSel / NumCountsSel gives the
// bin-averaged <gamma_t^1h, full>(R) ready to be summed with the
// gamma_t^prj projection observable in the likelihood.
#include "models/n_operator_sel_t.hh"
#include "modules/num_counts_sel/lensing_weights.hh"
#include "utils/module_macros.hh"

using Shear1hMisSel =
  y3_cluster::NOperatorSelRadial<y3_cluster_sel_weights::Shear1hMisWeight>;

DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE(Shear1hMisSel)
