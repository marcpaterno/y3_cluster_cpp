#ifndef Y3CLUSTER_FITS_LOADER_HH
#define Y3CLUSTER_FITS_LOADER_HH

#include "interp_1d.hh"

#include <string>
#include <utility>
#include <vector>

namespace y3_cluster {
  using doubles = std::vector<double>;

  // Returns a pair of {zgrid, {pz1, pz2, pz3, pz4}}
  std::pair<doubles, std::vector<doubles>> load_pzsource_data(
    const std::string& filename);

  std::vector<Interp1D> load_pzsources(const std::string& filename);
} // namespace y3_cluster

#endif
