#ifndef Y3_CLUSTER_AVERAGE_SCI_HH
#define Y3_CLUSTER_AVERAGE_SCI_HH

#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "models/ez.hh"
#include "utils/primitives.hh"
#include "utils/interp_1d.hh"
#include "utils/datablock_reader.hh"

#include <memory>
#include <cmath>

namespace y3_cluster
{
  class AVERAGE_SCI_t {
  private:
    std::shared_ptr<Interp1D const> _sci;

  public:
    AVERAGE_SCI_t(std::shared_ptr<Interp1D const> sci): _sci(sci){}

    using doubles = std::vector<double>;

    explicit AVERAGE_SCI_t(cosmosis::DataBlock& sample)
      : _sci(std::make_shared<Interp1D const>(
          get_datablock<doubles>(sample, "average_sigma_crit_inv", "zlense"),
          get_datablock<doubles>(sample, "average_sigma_crit_inv", "sci_average")))
    {}

    double
    operator()(double zt) const
    { 
      return _sci->eval(zt);
    }
  };
}

#endif
