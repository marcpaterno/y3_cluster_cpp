// N_i[gamma_t^1h(R)] — 1-halo tangential shear, DeltaSigma_nfw * Sigma_crit^-1.
// Ratio Shear1hSel / NumCountsSel gives the bin-averaged <gamma_t^1h(R)>.
#include "models/n_operator_sel_t.hh"
#include "modules/num_counts_sel/lensing_weights.hh"
#include "utils/module_macros.hh"

using Shear1hSel =
  y3_cluster::NOperatorSelRadial<y3_cluster_sel_weights::Shear1hWeight>;

DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE(Shear1hSel)
