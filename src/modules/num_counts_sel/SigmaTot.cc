// N_i[Sigma_tot(R)] — total (1h + b*2h) projected surface density.
// Ratio SigmaTotSel / NumCountsSel gives the bin-averaged <Sigma_tot(R)>.
#include "models/n_operator_sel_t.hh"
#include "modules/num_counts_sel/lensing_weights.hh"
#include "utils/module_macros.hh"

using SigmaTotSel =
  y3_cluster::NOperatorSelRadial<y3_cluster_sel_weights::SigmaTotWeight>;

DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE(SigmaTotSel)
