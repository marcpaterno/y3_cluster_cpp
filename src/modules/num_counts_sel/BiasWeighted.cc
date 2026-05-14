// N_i[b] — halo-bias-weighted number counts. Ratio N_i[b]/N_i[1] gives b_eff.
//
// b(M, z) is read from the "haloModel" section (populated by
// y3_buzzard/haloModelCosmosis.py) as a 2-D lookup table.
#include "models/n_operator_sel_t.hh"
#include "utils/interp_2d.hh"
#include "utils/make_interp_2d.hh"
#include "utils/module_macros.hh"

#include "cosmosis/datablock/datablock.hh"

#include <optional>

namespace {
  struct WeightBias {
    static char const* module_label() { return "BiasWeightedSel"; }
    explicit WeightBias(cosmosis::DataBlock&) {}
    void
    set_sample(cosmosis::DataBlock& sample)
    {
      bias_.emplace(
        y3_cluster::make_Interp2D(sample, "haloModel", "lnM", "z", "bias"));
    }
    double
    operator()(double lnM, double zt) const
    {
      return bias_->clamp(lnM, zt);
    }
  private:
    std::optional<y3_cluster::Interp2D> bias_;
  };
} // namespace

using BiasWeightedSel = y3_cluster::NOperatorSelScalar<WeightBias>;

DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE(BiasWeightedSel)
