#ifndef Y3_CLUSTER_PZSOURCE_HH
#define Y3_CLUSTER_PZSOURCE_HH

#include <datablock_reader.hh>
#include <cmath>
#include "primitives.hh"

namespace y3_cluster {
  namespace {
    std::vector<std::shared_ptr<y3_cluster::Interp1D const>>
    get_pzsource_distrs(cosmosis::DataBlock& db)
    {
      std::vector<std::shared_ptr<y3_cluster::Interp1D const>> output{};

      std::vector<double> zs = get_datablock<std::vector<double>>(db, "cluster_abundance", "zsource");
      std::vector<double> pzs1 = get_datablock<std::vector<double>>(db, "cluster_abundance", "pzsource_bin1"),
                          pzs2 = get_datablock<std::vector<double>>(db, "cluster_abundance", "pzsource_bin2"),
                          pzs3 = get_datablock<std::vector<double>>(db, "cluster_abundance", "pzsource_bin3"),
                          pzs4 = get_datablock<std::vector<double>>(db, "cluster_abundance", "pzsource_bin4");

      output.push_back(std::make_shared<y3_cluster::Interp1D const>(zs, pzs1));
      output.push_back(std::make_shared<y3_cluster::Interp1D const>(zs, pzs2));
      output.push_back(std::make_shared<y3_cluster::Interp1D const>(zs, pzs3));
      output.push_back(std::make_shared<y3_cluster::Interp1D const>(zs, pzs4));

      return output;
    }
  }

  class PZSOURCE_t {
  public:
    // TODO take in pzsource distributions
    explicit PZSOURCE_t(std::vector<std::shared_ptr<y3_cluster::Interp1D const>> pzsources)
        : pzsources(pzsources)
    {}

    explicit PZSOURCE_t(cosmosis::DataBlock& sample)
      : pzsources(get_pzsource_distrs(sample))
    {}

    double
    operator()(std::size_t bin_no, double zs) const
    {
      return pzsources[bin_no]->eval(zs);
    }

    std::size_t
    nsources() const
    {
      return pzsources.size();
    }

  private:
    // TODO move this to a different pzsource_t implementation?
    std::vector<std::shared_ptr<y3_cluster::Interp1D const>> pzsources;
  };
}

#endif
