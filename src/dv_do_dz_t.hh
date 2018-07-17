#ifndef Y3_CLUSTER_DV_DO_DZ_T_HH
#define Y3_CLUSTER_DV_DO_DZ_T_HH

#include <datablock_reader.hh>
#include <read_vector.hh>
#include "ez.hh"
#include "interp_1d.hh"

#include <memory>
#include <vector>

namespace y3_cluster {
  class DV_DO_DZ_t {
  public:
    DV_DO_DZ_t(std::shared_ptr<Interp1D const> da, y3_cluster::EZ ezt, double h)
      : _da(da), _ezt(ezt), _h(h)
    {}

    using doubles = std::vector<double>;

    explicit DV_DO_DZ_t(cosmosis::DataBlock& sample)
      : _da([sample] () {
	    // FIXME: This should be fully done from cosmosis
	    // when array is understood!
	    // Currently taken from trivial_gamma_t.cc
	    std::cout << "Starting DV_DO_DZ construction" << std::endl;
	    auto identity = [](double x) { return x; };
	    auto const da_arr = read_vector("d_a.txt", identity);
	    auto const zz_da = read_vector("z_da.txt", identity);
	    std::cout << "Beginning da construction" << std::endl;
	    auto da = std::make_shared<Interp1D const>(zz_da, da_arr);
	    std::cout << "Successfully made da" << std::endl;
	    return da;
	    }())
      //: _da(std::make_shared<Interp1D const>(
      //    get_datablock<doubles>(sample, "DV_D0_DZ_params", "xs"),
      //    get_datablock<doubles>(sample, "DV_D0_DZ_params", "ys")))
      , _ezt(y3_cluster::EZ(get_datablock<double>(sample, "cosmological_parameters", "omega_m"),
                            get_datablock<double>(sample, "cosmological_parameters", "omega_l"),
                            get_datablock<double>(sample, "cosmological_parameters", "omega_k")))
      // FIXME: Check which h value this should be!
      , _h(get_datablock<double>(sample, "cosmological_parameters", "h0"))
      {
	    std::cout << "Finsihed DV_DO_DZ construction" << std::endl;
      }

    double
    operator()(double zt) const
    {
      // TODO: Check this, as it looks incorrect
      double const da_z = _da->eval(zt); // da_z needs to be in Mpc
      return 3000.0 * (1.0 + zt) * (1.0 + zt) * da_z*_h * da_z*_h / _ezt(zt);
    }

  private:
    std::shared_ptr<Interp1D const> _da;
    y3_cluster::EZ _ezt;
    double _h;
  };
}

#endif
