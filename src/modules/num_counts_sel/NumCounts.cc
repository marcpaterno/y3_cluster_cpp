// N_i[1] — number counts in (lambda_ob, z_ob) bin using the closed-form
// Costanzi-Y3 richness selection function. See src/models/n_operator_sel_t.hh.
#include "models/n_operator_sel_t.hh"
#include "utils/module_macros.hh"

#include "cosmosis/datablock/datablock.hh"

namespace {
  struct WeightOne {
    static char const* module_label() { return "NumCountsSel"; }
    explicit WeightOne(cosmosis::DataBlock&) {}
    void set_sample(cosmosis::DataBlock&) {}
    double operator()(double /*lnM*/, double /*zt*/) const { return 1.0; }
  };
} // namespace

using NumCountsSel = y3_cluster::NOperatorSelScalar<WeightOne>;

DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE(NumCountsSel)
