// Fixed-GL driver for the three Costanzi-2026 P[X] scalars.  One
// `set_sample` pass co-computes (P1, I1, J) on the common (zob_low,
// zob_high, lambda_bin) wall grid with J = I2 - I1 computed directly
// as int b (1 - sigma(theta)) xi f_A ... dz dlnM dtheta.  Computing
// J directly avoids the catastrophic cancellation at large theta
// where sigma(theta) -> 1 and I2 - I1 suffers precision loss; J is
// the denominator of b_zero in bias_from_precomp and is guaranteed
// >= 0 when xi > 0.  The evaluator writes three separate `vals`
// arrays to sections b_sel_marg_P1 / _I1 / _J consumed by bsel.py.
#include "models/p_operator_t.hh"
#include "utils/module_macros.hh"

using BSelMargIntegrand = y3_cluster::P_operator;

DEFINE_COSMOSIS_SCALAR_EVALUATOR_MODULE(BSelMargIntegrand)
