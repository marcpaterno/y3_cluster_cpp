#ifndef Y3_CLUSTER_Y_SZ
#define Y3_CLUSTER_Y_SZ

#include "cosmosis/datablock/ndarray.hh"
#include "utils/datablock_reader.hh"
#include "utils/interp_2d.hh"

#include <memory>

namespace y3_cluster {
  class Y_SZ {
    // TODO add members
    std::vector<std::shared_ptr<y3_cluster::Interp2D>> interp_1h;
    std::vector<std::shared_ptr<y3_cluster::Interp2D>> interp_2h;

    static bool
    valid_ys_arr(const std::vector<double>& ms,
                 const std::vector<double>& zs,
                 const std::vector<double>& thetas,
                 const cosmosis::ndarray<double>& ys)
    {
      auto extents = ys.extents();
      if (extents.size() != 3)
        return false;
      return (extents[0] == thetas.size()) & (extents[1] == ms.size()) &
             (extents[2] == zs.size());
    }

    static std::shared_ptr<y3_cluster::Interp2D>
    make_interp(const std::vector<double>& ms,
                const std::vector<double>& zs,
                // No const because the compiler gets confused
                // but this should be const
                cosmosis::ndarray<double>& ys_3d,
                unsigned theta_i)
    {
      std::vector<double> ys;
      for (auto i = 0u; i < zs.size(); i++) {
        for (auto j = 0u; j < ms.size(); j++) {
          ys.push_back(ys_3d(theta_i, j, i));
        }
      }
      return std::make_shared<y3_cluster::Interp2D>(ms, zs, ys);
    }

  public:
    explicit Y_SZ(cosmosis::DataBlock& sample) : interp_1h(), interp_2h()
    {
      using doubles = std::vector<double>;
      auto ms = get_datablock<doubles>(sample, "cluster_sz", "ms");
      auto zs = get_datablock<doubles>(sample, "cluster_sz", "zs");
      auto thetas = get_datablock<doubles>(sample, "cluster_sz", "thetas");

      auto ys_1h =
        get_datablock<cosmosis::ndarray<double>>(sample, "cluster_sz", "ys_1h");
      auto ys_2h =
        get_datablock<cosmosis::ndarray<double>>(sample, "cluster_sz", "ys_2h");

      if (!valid_ys_arr(ms, zs, thetas, ys_1h))
        throw std::invalid_argument("Invalid 1halo dimensions");
      if (!valid_ys_arr(ms, zs, thetas, ys_2h))
        throw std::invalid_argument("Invalid 2halo dimensions");

      for (auto ti = 0u; ti < thetas.size(); ti++) {
        interp_1h.push_back(make_interp(ms, zs, ys_1h, ti));
        interp_2h.push_back(make_interp(ms, zs, ys_2h, ti));
      }

      if (interp_1h.size() != interp_2h.size())
        throw std::runtime_error("Internal error: 1h and 2h terms"
                                 "have different dimensions");
    }

    // In place of attempting to do 3D interpolation, we operate on
    // discrete, zero-indexed angular bins
    double
    operator()(double M, double z, unsigned theta_bin) const
    {
      return interp_1h[theta_bin]->eval(M, z) +
             interp_2h[theta_bin]->eval(M, z);
    }
  };
}

#endif
