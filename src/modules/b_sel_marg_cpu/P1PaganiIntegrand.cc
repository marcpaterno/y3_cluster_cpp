// P1 operator integrand (Costanzi-2026 eq. 4, F = 1) -- adaptive Cuhre
// variant.  Kept for benchmarking vs the fixed-GL evaluator in
// P1Integrand.cc.
#include "models/p_operator_cuhre_t.hh"
#include "utils/module_macros.hh"

using BSelMargP1PaganiIntegrand =
    y3_cluster::P_operator_cuhre<y3_cluster::BSelMargP1PaganiWeight>;

DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE(BSelMargP1PaganiIntegrand)
