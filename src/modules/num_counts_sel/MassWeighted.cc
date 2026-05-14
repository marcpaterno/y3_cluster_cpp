// N_i[M] — mass-weighted number counts. Ratio N_i[M]/N_i[1] gives M_eff.
#include "models/n_operator_sel_t.hh"
#include "utils/module_macros.hh"

#include "cosmosis/datablock/datablock.hh"

#include <cmath>

namespace {
  struct WeightMass {
    static char const* module_label() { return "MassWeightedSel"; }
    explicit WeightMass(cosmosis::DataBlock&) {}
    void set_sample(cosmosis::DataBlock&) {}
    double operator()(double lnM, double /*zt*/) const { return std::exp(lnM); }
  };
} // namespace

using MassWeightedSel = y3_cluster::NOperatorSelScalar<WeightMass>;

DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE(MassWeightedSel)
