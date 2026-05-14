// N_i[DeltaSigma_1h(R)] — 1-halo differential surface density, R-by-R.
// Ratio DSigma1hSel / NumCountsSel gives the bin-averaged <DeltaSigma_1h(R)>.
#include "models/n_operator_sel_t.hh"
#include "modules/num_counts_sel/lensing_weights.hh"
#include "utils/module_macros.hh"

using DSigma1hSel =
  y3_cluster::NOperatorSelRadial<y3_cluster_sel_weights::DSigma1hWeight>;

DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE(DSigma1hSel)
