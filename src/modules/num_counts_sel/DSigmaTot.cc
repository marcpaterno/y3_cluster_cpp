// N_i[DeltaSigma_tot(R)] — total (1h + b*2h) differential surface density.
// Ratio DSigmaTotSel / NumCountsSel gives the bin-averaged <DeltaSigma_tot(R)>.
#include "models/n_operator_sel_t.hh"
#include "modules/num_counts_sel/lensing_weights.hh"
#include "utils/module_macros.hh"

using DSigmaTotSel =
  y3_cluster::NOperatorSelRadial<y3_cluster_sel_weights::DSigmaTotWeight>;

DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE(DSigmaTotSel)
