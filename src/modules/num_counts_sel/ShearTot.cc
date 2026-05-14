// N_i[gamma_t^tot(R)] — total-halo tangential shear, DeltaSigma_tot * Sigma_crit^-1.
// Ratio ShearTotSel / NumCountsSel gives the bin-averaged <gamma_t^tot(R)>.
#include "models/n_operator_sel_t.hh"
#include "modules/num_counts_sel/lensing_weights.hh"
#include "utils/module_macros.hh"

using ShearTotSel =
  y3_cluster::NOperatorSelRadial<y3_cluster_sel_weights::ShearTotWeight>;

DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE(ShearTotSel)
