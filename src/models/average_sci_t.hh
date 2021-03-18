#ifndef Y3_CLUSTER_AVERAGE_SCI_HH
#define Y3_CLUSTER_AVERAGE_SCI_HH

#include "cosmosis/datablock/datablock.hh"
#include "utils/interp_1d.hh"

namespace y3_cluster {
  class AVERAGE_SCI_t {
  private:
    Interp1D _sci;

  public:
    AVERAGE_SCI_t(Interp1D const& sci) : _sci(sci) {}

    using doubles = std::vector<double>;

    explicit AVERAGE_SCI_t(cosmosis::DataBlock& config)
      : _sci(config.view<doubles>("average_sigma_crit_inv", "zlense"),
             config.view<doubles>("average_sigma_crit_inv", "sci_average"))
    {}

    double
    operator()(double zt) const
    {
      return _sci(zt);
    }
  };
}

#endif
