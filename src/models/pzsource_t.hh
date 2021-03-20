#ifndef Y3_CLUSTER_PZSOURCE_HH
#define Y3_CLUSTER_PZSOURCE_HH

#include <stdexcept>

#include "cosmosis/datablock/datablock.hh"
#include "utils/interp_1d.hh"
#include "utils/make_interp_1d.hh"
#include "utils/primitives.hh"

namespace y3_cluster {
  class PZSOURCE_t {
  public:

    // TODO take in pzsource distributions
    explicit PZSOURCE_t(
      std::vector<y3_cluster::Interp1D> const& pzsources,
      std::vector<double> const&  zoffsets)
      : pzsources(pzsources), zoffsets(zoffsets)
    {}

    explicit PZSOURCE_t(
      cosmosis::DataBlock& sample, std::vector<y3_cluster::Interp1D> pzsources)
      : pzsources(pzsources)
      , zoffsets{sample.view<double>("cluster_abundance", "pzsource_offset_1"),
                 sample.view<double>("cluster_abundance", "pzsource_offset_2"),
                 sample.view<double>("cluster_abundance", "pzsource_offset_3"),
                 sample.view<double>("cluster_abundance", "pzsource_offset_4")}
    {
      if (pzsources.size() != zoffsets.size())
        throw std::runtime_error(
          "number of zsource distributions and number of offsets do not match");
    }

    double
    operator()(std::size_t bin_no, double zs) const
    {
      return pzsources[bin_no](zs + zoffsets[bin_no]);
    }

    std::size_t
    nsources() const
    {
      return pzsources.size();
    }

  private:
    std::vector<y3_cluster::Interp1D> pzsources;
    std::vector<double> zoffsets;
  };
}

#endif
