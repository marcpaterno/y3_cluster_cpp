#include "cosmosis/datablock/ndarray.hh"
#include "utils/datablock_reader.hh"
#include "utils/make_grid_points.hh"
#include "utils/make_integration_volumes.hh"
#include "utils/module_macros.hh"

#include "cosmosis/datablock/datablock.hh"
#include "cubacpp/integration_result.hh"
#include "cubacpp/integration_volume.hh"

#include "models/hmf_t.hh"
#include "models/lc_lt_y1_t.hh"
#include "models/sptxdes/mor_y3xspt_t.hh"

#include <cmath>
#include <optional>
#include <string>
#include <vector>

using cosmosis::DataBlock;

class AbundanceIntegralDESxSPT {
public:
  // 4 varied grid parameters: lo, xi, zo, gamma_field
  using grid_t = y3_cluster::grid_t<4>;
  using grid_point_t = grid_t::value_type;

private:
  // Number of integration dimensions
  using volume_t = cubacpp::IntegrationVolume<2>;

  // State obtained from configuration
  y3_cluster::LC_LT_Y1_t PLamobs_LamtrueZtrue;

  // State obtained from each sample.
  // If there were a type X that did not have a default constructor,
  // we would use std::optional<X> as our data member.
  std::optional<y3_cluster::HMF_t> hmf;
  std::optional<y3_cluster::MOR_Y3xSPT_t> PLamtrueZeta_lnMZtrue;

  // State set for current 'bin' to be integrated.
  double lamobs_;
  double zeta_;
  double zobs_;
  double gamma_field_;

public:
  // Initialize my integrand object from the parameters read
  // from the relevant block in the CosmoSIS ini file.
  explicit AbundanceIntegralDESxSPT(DataBlock&);

  // Set any data members from values read from the current sample.
  // Do not attempt to copy the sample!.
  void set_sample(DataBlock& sample);

  // Set the data for the current bin.
  void set_grid_point(grid_point_t const& grid_point);

  // The function to be integrated. All arguments to this function must be of
  // type double, and there must be at least two of them (because our
  // integration routine does not work for functions of one variable). The
  // function is const because calling it does not change the state of the
  // object.
  double operator()(double lamtrue, double lnM200m) const;

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
  static std::vector<volume_t> make_integration_volumes(DataBlock& cfg);

  // The following non-member (static) function creates a vector of grid points
  // on which the integration results are to be evaluated, based on parameters
  // read from the configuration block for the module.
  static grid_t make_grid_points(DataBlock& cfg);
};

// Set up integral that does not depend on sample
AbundanceIntegralDESxSPT::AbundanceIntegralDESxSPT(DataBlock&)
  : PLamobs_LamtrueZtrue()
  , hmf()
  , PLamtrueZeta_lnMZtrue()
  , lamobs_()
  , zeta_()
  , zobs_()
  , gamma_field_()
{}

// Set up the integration using the current sample
void
AbundanceIntegralDESxSPT::set_sample(DataBlock& sample)
{
  hmf.emplace(sample);
  PLamtrueZeta_lnMZtrue.emplace(sample);
}

// Set the grid points for the current integration
void
AbundanceIntegralDESxSPT::set_grid_point(grid_point_t const& grid_point)
{
  lamobs_ = grid_point[0];
  zeta_ = grid_point[1];
  zobs_ = grid_point[2];
  gamma_field_ = grid_point[3];
}

// Evaluate the integrand
// Anything initilized with std::optional needs an asterisk
double
AbundanceIntegralDESxSPT::operator()(double lamtrue, double lnM200m) const
{
  return (*hmf)(lnM200m, zobs_) *
         PLamobs_LamtrueZtrue(lamobs_, lamtrue, zobs_) *
         (*PLamtrueZeta_lnMZtrue)(lamtrue, zeta_, zobs_, lnM200m, gamma_field_);
}

// Name of the module in the datablock
char const*
AbundanceIntegralDESxSPT::module_label()
{
  return "abundance_integral";
}

// Set the bounds of integration
std::vector<AbundanceIntegralDESxSPT::volume_t>
AbundanceIntegralDESxSPT::make_integration_volumes(cosmosis::DataBlock& cfg)
{
  // Make some handy definitions
  auto modulelabel = AbundanceIntegralDESxSPT::module_label();
  std::vector<std::string> dimensions{"lamtrue", "lnM200m"};
  std::size_t const nvolumes = cfg.view<int>(modulelabel, "nclusters");
  std::size_t const N = 2;

  // Initialize the arrays of bounds
  std::vector<cubacpp::array<N>> lowbounds(nvolumes);
  std::vector<cubacpp::array<N>> highbounds(nvolumes);

  // Load lamtrue bounds and set to first integration dimension
  auto lamtrue_bounds = get_vector_double(cfg, modulelabel, "lamtrue_bounds");
  for (std::size_t ivol = 0; ivol != nvolumes; ++ivol) {
    lowbounds[ivol][0] = lamtrue_bounds[0];
    highbounds[ivol][0] = lamtrue_bounds[1];
  }

  // Make lnM200m the second integration dimension
  auto lnM200m_bounds = get_vector_double(cfg, modulelabel, "lnM200m_bounds");
  for (std::size_t ivol = 0; ivol != nvolumes; ++ivol) {
    lowbounds[ivol][1] = lnM200m_bounds[0];
    highbounds[ivol][1] = lnM200m_bounds[1];
  }

  // Convert into the expected IntegrationVolume object
  std::vector<cubacpp::IntegrationVolume<N>> result;
  result.reserve(lowbounds.size());
  for (std::size_t i = 0; i != lowbounds.size(); ++i) {
    result.emplace_back(lowbounds[i], highbounds[i]);
  }
  return result;
}

// Set the grid points to evaluate the integral at
AbundanceIntegralDESxSPT::grid_t
AbundanceIntegralDESxSPT::make_grid_points(cosmosis::DataBlock& cfg)
{
  return y3_cluster::load_grid_from_file_wall_of_numbers(
    cfg,
    AbundanceIntegralDESxSPT::module_label(),
    "lamobs",
    "xi",
    "zobs",
    "gamma_field");
}

DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE(AbundanceIntegralDESxSPT)
