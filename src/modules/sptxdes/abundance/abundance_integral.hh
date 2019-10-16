#ifndef Y3_CLUSTER_ABUNDANCE_INTEGRAL_H
#define Y3_CLUSTER_ABUNDANCE_INTEGRAL_H

#include "utils/make_integration_volumes.hh"
#include "utils/module_macros.hh"

#include "cosmosis/datablock/datablock.hh"
#include "cubacpp/integration_result.hh"
#include "cubacpp/integration_volume.hh"

#include "models/zo_zt_des_t.hh"
#include "models/omega_z_y3xspt.hh"
#include "models/dv_do_dz_t.hh"
#include "models/hmf_t.hh"
#include "models/lc_lt_t.hh"
#include "models/pxizeta_t.hh"
#include "models/mor_y3xspt_t.hh"

#include <optional>
#include <vector>

using namespace y3_cluster;

class AbundanceIntegralDesXSpt{
public:
  // 3 varied grid parameters: lo, xi, zo
  using grid_point_t = std::array<double, 3>;

private:
  // Number of integration dimensions
  using volume_t = cubacpp::IntegrationVolume<4>;

  // State obtained from configuration
  ZO_ZT_DES_t p_zo_zt;
  OMEGA_Z_Y3XSPT omega_z;
  LC_LT_t p_lo_ltzt;
  PXiZeta_t p_xi_zeta;

  // State obtained from each sample.
  // If there were a type X that did not have a default constructor,
  // we would use std::optional<X> as our data member.
  std::optional<DV_DO_DZ_t> dv_dodz;
  std::optional<HMF_t> hmf;
  std::optional<MOR_Y3xSPT_t> p_ltzeta_lnmzt;

  // State set for current 'bin' to be integrated.
  double lamobs_;
  double xi_;
  double zobs_;
  double gamma_field_;

public:
  // Initialize my integrand object from the parameters read
  // from the relevant block in the CosmoSIS ini file.
  explicit AbundanceIntegralDesXSpt(cosmosis::DataBlock& config);

  // Set any data members from values read from the current sample.
  // Do not attempt to copy the sample!.
  void set_sample(cosmosis::DataBlock& sample);

  // Set the data for the current bin.
  void set_grid_point(grid_point_t const& grid_point);

  // The function to be integrated. All arguments to this function must be of
  // type double, and there must be at least two of them (because our
  // integration routine does not work for functions of one variable). The
  // function is const because calling it does not change the state of the
  // object.
  double operator()(double ztrue, double lnM200m, double lamtrue, double zeta) const;

  // module_label() is a non-member (static) function that returns the label for
  // this module. The name this returns
  // is the name that must be used in the 'ini file' for configuring the module
  // made with this class.
  // We return char const* rather than std::string to avoid some needless memory
  // allocations.
  static char const* module_label();

  // The following non-member (static) function creates a vector of integration
  // volumes (the type alias defined above) based on the parameters read from
  // the configuration block for the module.
  static std::vector<volume_t> make_integration_volumes(
    cosmosis::DataBlock& cfg);

  // The following non-member (static) function creates a vector of grid points
  // on which the integration results are to be evaluated, based on parameters
  // read from the configuration block for the module.
  static std::vector<grid_point_t> make_grid_points(cosmosis::DataBlock& cfg);
};

#endif
