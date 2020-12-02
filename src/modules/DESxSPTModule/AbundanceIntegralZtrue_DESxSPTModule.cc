#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/ndarray.hh"
#include "utils/datablock_reader.hh"
#include "utils/make_grid_points.hh"
#include "utils/make_integration_volumes.hh"
#include "utils/module_macros.hh"
#include "utils/OneDIntegrationModule.hh"

#include "models/dv_do_dz_t.hh"
#include "models/sptxdes/omega_z_y3xspt.hh"
#include "models/zo_zt_des_t.hh"

#include <cmath>
#include <optional>
#include <string>
#include <vector>

using cosmosis::DataBlock;

class ZtrueIntegralDESxSPT {
public:
  // The number of grid dimensions. I only use zobs but
  // load all of the cluster observables to have a single
  // catalog file without making the loading code complicated
  using grid_t = y3_cluster::grid_t<4>;
  using grid_point_t = grid_t::value_type;

private:
  // From configuration - constructor
  y3_cluster::ZO_ZT_DES_t PZobsZtrue;
  y3_cluster::OMEGA_Z_Y3XSPT OmegaZtrue;

  // From likelihood sample - sample
  std::optional<y3_cluster::DV_DO_DZ_t> dv_dodzt;

  // From cluster sample - grid
  double zobs_;

public:
  // Module name
  static char const* module_label();

  // Integration variable name
  static char const* integration_variable();

  // The integrand
  double operator()(double ztrue) const;

  // Initialize from the config file
  explicit ZtrueIntegralDESxSPT(DataBlock&);

  // Update the object using the current likelihood sample
  void set_sample(DataBlock& sample);

  // Set the data for the current bin.
  void set_grid_point(grid_point_t const& pt);

  // The following non-member (static) function creates a vector of grid points
  // on which the integration results are to be evaluated, based on parameters
  // read from the configuration block for the module.
  static grid_t make_grid_points(DataBlock& cfg);

  // The following non-member (static) function creates a vector of integration
  // volumes (the type alias defined above) based on the parameters read from
  // the configuration block for the module.
  static std::vector<std::array<double, 2>> make_integration_volumes(
    DataBlock& cfg);
};

ZtrueIntegralDESxSPT::ZtrueIntegralDESxSPT(DataBlock&)
  : PZobsZtrue(), OmegaZtrue(), dv_dodzt(), zobs_()
{}

char const*
ZtrueIntegralDESxSPT::module_label()
{
  return "ztrue_integration";
}

char const*
ZtrueIntegralDESxSPT::integration_variable()
{
  return "ztrue";
}

double
ZtrueIntegralDESxSPT::operator()(double ztrue) const
{
  return PZobsZtrue(zobs_, ztrue) * OmegaZtrue(ztrue) * (*dv_dodzt)(ztrue);
}

void
ZtrueIntegralDESxSPT::set_sample(DataBlock& sample)
{
  dv_dodzt.emplace(sample);
}

void
ZtrueIntegralDESxSPT::set_grid_point(grid_point_t const& grid_point)
{
  zobs_ = grid_point[2];
}

ZtrueIntegralDESxSPT::grid_t
ZtrueIntegralDESxSPT::make_grid_points(DataBlock& cfg)
{
  return y3_cluster::load_grid_from_file_wall_of_numbers(
    cfg, module_label(), "lamobs", "xi", "zobs", "gamma_field");
}

std::vector<std::array<double, 2>>
ZtrueIntegralDESxSPT::make_integration_volumes(DataBlock& cfg)
{
  // Load zobs from gridfile
  auto modulelabel = module_label();
  auto gridpts = y3_cluster::load_grid_from_file_wall_of_numbers(
    cfg, modulelabel, "lamobs", "xi", "zobs", "gamma_field");
  std::vector<double> zobs;
  for (std::size_t i = 0; i < gridpts.size(); ++i)
    zobs.push_back(gridpts[i][2]);

  // Load from the config file
  std::size_t const nvolumes = cfg.view<int>(modulelabel, "nclusters");
  double Nsigma = cfg.view<double>(modulelabel, "Nsigma_ztrue");

  // Bounds on ztrue are zobs +- Nsigma_ztrue*sigma_z(zobs)
  std::vector<double> lows, highs;
  y3_cluster::SIGMA_PHOTOZ_DES_t sigma_model;
  for (std::size_t ivol = 0; ivol != nvolumes; ++ivol) {
    lows.push_back(zobs[ivol] - Nsigma * sigma_model(zobs[ivol]));
    highs.push_back(zobs[ivol] + Nsigma * sigma_model(zobs[ivol]));
  }

  // Put into an std:array and return
  std::vector<std::array<double, 2>> volumes;
  for (std::size_t i = 0; i != lows.size(); ++i)
    volumes.push_back(std::array<double, 2>{lows[i], highs[i]});
  return volumes;
}

DEFINE_COSMOSIS_ONED_INTEGRATION_MODULE(ZtrueIntegralDESxSPT)
