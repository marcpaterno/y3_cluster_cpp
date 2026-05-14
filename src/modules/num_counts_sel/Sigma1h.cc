// N_i[Sigma_1h(R)] — 1-halo projected surface density, R-by-R.
// Ratio Sigma1hSel / NumCountsSel gives the bin-averaged <Sigma_1h(R)>.
#include "models/n_operator_sel_t.hh"
#include "modules/num_counts_sel/lensing_weights.hh"
#include "utils/module_macros.hh"

using Sigma1hSel =
  y3_cluster::NOperatorSelRadial<y3_cluster_sel_weights::Sigma1hWeight>;

DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE(Sigma1hSel)
