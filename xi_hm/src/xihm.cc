#include "/cosmosis/cosmosis/datablock/datablock.hh"
#include <datablock_reader.hh>
#include <nfw.hh>
#include <xinl.hh>
#include <del_sig_y3.hh>

namespace y3_cluster {
  class DeltaSigmaModule {
    double ln_M_min;
    double ln_M_max;
    std::size_t M_bins;
    double z_min;
    double z_max;
    std::size_t z_bins;
    double ln_R_perp_min;
    double ln_R_perp_max;
    std::size_t R_perp_bins;

  public:
    explicit DeltaSigmaModule(cosmosis::DataBlock& config);
    DATABLOCK_STATUS execute(cosmosis::DataBlock& sample);
  };
}

using namespace y3_cluster;


y3_cluster::DeltaSigmaModule::DeltaSigmaModule(cosmosis::DataBlock& config)
{
	ln_M_min = std::log( get_datablock<double>(config, OPTION_SECTION, "M_min") ); 
	ln_M_max = std::log( get_datablock<double>(config, OPTION_SECTION, "M_max") ); 
	M_bins = get_datablock<int>(config, OPTION_SECTION, "M_bins"); 
	z_min = get_datablock<double>(config, OPTION_SECTION, "z_min") ; 
	z_max = get_datablock<double>(config, OPTION_SECTION, "z_max") ; 
	z_bins = get_datablock<int>(config, OPTION_SECTION, "z_bins"); 

	ln_R_perp_min = std::log( get_datablock<double>(config, OPTION_SECTION, "R_perp_min") ); 
	ln_R_perp_max = std::log( get_datablock<double>(config, OPTION_SECTION, "R_perp_max") ); 
	R_perp_bins = get_datablock<int>(config, OPTION_SECTION, "R_perp_bins"); 

}


DATABLOCK_STATUS
y3_cluster::DeltaSigmaModule::execute(cosmosis::DataBlock& sample)
{

  // Pull arrays from the data block:
  const double concentration = 5.0;
  nfw_dsigma nfw(concentration);
  xi_nl xinl(sample);
  std::vector<double>  xi_1halo(R_perp_bins * M_bins);
  std::vector<double>  xi_2halo(R_perp_bins * z_bins);

  // set up return arrays
  auto R_perp = std::vector<double> (R_perp_bins);
  for (std::size_t i = 0; i < R_perp_bins; i++)
	      R_perp[i] = std::exp( (ln_R_perp_max - ln_R_perp_min ) / R_perp_bins * i + ln_R_perp_min );

  auto lnM = std::vector<double> (M_bins);
  auto m_h = std::vector<double> (M_bins);
  auto z = std::vector<double> (z_bins);
  for (std::size_t i = 0; i < M_bins; i++) {
	      lnM[i] = (ln_M_max - ln_M_min) / M_bins * i + ln_M_min;
	      m_h[i] = std::exp(lnM[i]);
   }
  for (std::size_t i = 0; i < z_bins; i++) 
	      z[i] = (z_max - z_min) / z_bins * i + z_min;
  for (std::size_t i = 0; i < R_perp_bins; i++) {
	  for (std::size_t j = 0; j < M_bins; j++) {
		xi_1halo[i*M_bins+j] = nfw(R_perp[i], m_h[j]);
	  }
	  for (std::size_t kk = 0; kk < z_bins; kk++) {
		xi_2halo[i*z_bins+kk] = xinl(R_perp[i], z[kk]);
	  }
  }


  // Store arrays in the datablock
  // Need to add cosmosis datablock sections
  std::vector<std::size_t> extents_m{{R_perp_bins, M_bins}};
  std::vector<std::size_t> extents_z{{R_perp_bins, z_bins}};
  sample.put_val<cosmosis::ndarray<double>>("cluster_abundance", "xi_nfw", {xi_1halo, extents_m});
  sample.put_val<cosmosis::ndarray<double>>("cluster_abundance", "xi_nl", {xi_2halo, extents_z});
  sample.put_val<std::vector<double>>("cluster_abundance", "R_perp", R_perp);
  sample.put_val<std::vector<double>>("cluster_abundance", "m_h", m_h);
  sample.put_val<std::vector<double>>("cluster_abundance", "z", z);
  sample.put_val<std::vector<double>>("cluster_abundance", "lnM", lnM);
  
  DEL_SIG_y3 ds(sample);
  double temp_m, res;
  for (std::size_t i = 0; i < 500; i++) {
    temp_m = (33.0 - 32.0) / 500.0 * i + 32.0;
    res= ds(0.3, temp_m, 0.3);
				         }
  return DBS_SUCCESS;
}



// Module interface implementation

extern "C" {

	void *
	setup(cosmosis::DataBlock *options)
	{
	    return new DeltaSigmaModule(*options);
	}

	DATABLOCK_STATUS
	execute(cosmosis::DataBlock *block, void *config)
	{
	    auto mod = static_cast<DeltaSigmaModule*>(config);
	    return mod->execute(*block);
	}

	int
	cleanup(void *config)
	{
	    delete static_cast<DeltaSigmaModule*>(config);
	    return 0;
	}

}
