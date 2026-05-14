// I2 operator integrand -- adaptive Cuhre variant (benchmark).
#include "models/p_operator_cuhre_t.hh"
#include "utils/module_macros.hh"

using BSelMargI2PaganiIntegrand =
    y3_cluster::P_operator_cuhre<y3_cluster::BSelMargI2PaganiWeight>;

DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE(BSelMargI2PaganiIntegrand)
