#include "cosmosis/datablock/datablock.hh"
#include "utils/datablock_reader.hh"
#include "models/xinl.hh"
#include "models/xi_to_sigma_t.hh"
#include "models/nfw_xi_t.hh"
#include "models/nfw_sigma_t.hh"

namespace y3_cluster {
  class SigmaModule {
    double ln_M_min;
    double ln_M_max;
    std::size_t M_bins;
    double z_min;
    double z_max;
    std::size_t z_bins;
    double ln_R2d_min;
    double ln_R2d_max;
    std::size_t R2d_bins;
    double ln_R3d_min;
    double ln_R3d_max;
    std::size_t R3d_bins;

  public:
    explicit SigmaModule(cosmosis::DataBlock& config);
    DATABLOCK_STATUS execute(cosmosis::DataBlock& sample);
  };
}

using namespace y3_cluster;


y3_cluster::SigmaModule::SigmaModule(cosmosis::DataBlock& config)
{
	ln_M_min = std::log( get_datablock<double>(config, OPTION_SECTION, "M_min") ); 
	ln_M_max = std::log( get_datablock<double>(config, OPTION_SECTION, "M_max") ); 
	M_bins = get_datablock<int>(config, OPTION_SECTION, "M_bins"); 
	z_min = get_datablock<double>(config, OPTION_SECTION, "z_min") ; 
	z_max = get_datablock<double>(config, OPTION_SECTION, "z_max") ; 
	z_bins = get_datablock<int>(config, OPTION_SECTION, "z_bins"); 

	ln_R2d_min = std::log( get_datablock<double>(config, OPTION_SECTION, "R_perp_min") ); 
	ln_R2d_max = std::log( get_datablock<double>(config, OPTION_SECTION, "R_perp_max") ); 
	R2d_bins = get_datablock<int>(config, OPTION_SECTION, "R_perp_bins"); 

	ln_R3d_min = std::log( get_datablock<double>(config, OPTION_SECTION, "Radii_min") ); 
	ln_R3d_max = std::log( get_datablock<double>(config, OPTION_SECTION, "Radii_max") ); 
	R3d_bins = get_datablock<int>(config, OPTION_SECTION, "Radii_bins"); 
}


DATABLOCK_STATUS
y3_cluster::SigmaModule::execute(cosmosis::DataBlock& sample)
{

  // set up classes:
  xi_nl xinl(sample);
  std::vector<double>  xi_2halo(R3d_bins * z_bins);
  std::vector<double>  xi_1halo(R3d_bins * M_bins);
  std::vector<double>  sigma_1halo(R3d_bins * M_bins);
  //

  // set up radii
  auto R2d = std::vector<double> (R2d_bins);
  for (std::size_t i = 0; i < R2d_bins; i++)
	      R2d[i] = std::exp( (ln_R2d_max - ln_R2d_min ) / R2d_bins * i + ln_R2d_min );
  auto R3d = std::vector<double> (R3d_bins);
  for (std::size_t i = 0; i < R3d_bins; i++)
	      R3d[i] = std::exp( (ln_R3d_max - ln_R3d_min ) / R3d_bins * i + ln_R3d_min );
  //

  auto z= std::vector<double> (z_bins);
  for (std::size_t i = 0; i < z_bins; i++) 
	      z[i] = (z_max - z_min) / z_bins * i + z_min;
  auto lnM = std::vector<double> (M_bins);
  auto m_h = std::vector<double> (M_bins);
  for (std::size_t i = 0; i < M_bins; i++) {
          lnM[i] = (ln_M_max - ln_M_min) / M_bins * i + ln_M_min;
          m_h[i] = std::exp(lnM[i]);
   }

  for (std::size_t i = 0; i < R3d_bins; i++) 
  {
	  for (std::size_t kk = 0; kk < z_bins; kk++) 
      {
		xi_2halo[i+kk*R3d_bins] = xinl(R3d[i], z[kk]);
	  }
  }
   
  nfw_xi nfw(sample);
  nfw_sigma nfw_sigma(sample);
  for (std::size_t i = 0; i < R3d_bins; i++) 
  {
	  for (std::size_t kk = 0; kk < M_bins; kk++) 
      {
		xi_1halo[i+kk*R3d_bins] = nfw(R3d[i], m_h[kk]);
	  }
  }
  for (std::size_t i = 0; i < R2d_bins; i++) 
  {
	  for (std::size_t kk = 0; kk < M_bins; kk++) 
      {
		sigma_1halo[i+kk*R3d_bins] = nfw_sigma(R2d[i], m_h[kk]);
	  }
  }


  // Store arrays in the datablock
  // Need to add cosmosis datablock sections
  std::vector<std::size_t> extents_3d{{z_bins, R3d_bins}};
  sample.put_val<cosmosis::ndarray<double>>("deltasigma", "Xi_2", {xi_2halo, extents_3d});
  sample.put_val<std::vector<double>>("deltasigma", "R_Xi", R3d);
  sample.put_val<std::vector<double>>("deltasigma", "z", z);

  std::vector<std::size_t> extents_3d_1{{M_bins, R3d_bins}};
  sample.put_val<cosmosis::ndarray<double>>("deltasigma", "Xi_1", {xi_1halo, extents_3d_1});
  sample.put_val<std::vector<double>>("deltasigma", "lnM", lnM);
  sample.put_val<std::vector<double>>("deltasigma", "m_h", m_h);
  std::vector<std::size_t> extents_2d_1{{M_bins, R2d_bins}};
  sample.put_val<cosmosis::ndarray<double>>("deltasigma", "sigma_1", {sigma_1halo, extents_2d_1});
  
  SIG_y3 sig(z, R3d, xi_2halo, sample);
  std::vector<double>  sigma_2halo(R2d_bins * z_bins);
  for (std::size_t i = 0; i < R2d_bins; i++)
  {
	  for (std::size_t kk = 0; kk < z_bins; kk++) 
      {
          sigma_2halo[i+kk*R2d_bins] = sig(R2d[i],z[kk]);
      }
  } 
  std::vector<std::size_t> extents_2d{{z_bins, R2d_bins}};
  sample.put_val<cosmosis::ndarray<double>>("deltasigma", "sigma_2", {sigma_2halo, extents_2d});
  sample.put_val<std::vector<double>>("deltasigma", "R_sigma_deltasigma", R2d);

  

  return DBS_SUCCESS;
}



// Module interface implementation

extern "C" {

	void *
	setup(cosmosis::DataBlock *options)
	{
	    return new SigmaModule(*options);
	}

	DATABLOCK_STATUS
	execute(cosmosis::DataBlock *block, void *config)
	{
	    auto mod = static_cast<SigmaModule*>(config);
	    return mod->execute(*block);
	}

	int
	cleanup(void *config)
	{
	    delete static_cast<SigmaModule*>(config);
	    return 0;
	}

}
