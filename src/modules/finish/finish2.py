import numpy as np
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.interpolate import interp1d

##The aim of this module is to convert sigma into deltasigma
#this is made with changing the shape
def setup(options):
   
    section = option_section
    
   
    # model module     
    rad = options[section, "radius"]
    
    
    #length of radii
    n_rad = len(rad)


    return rad, n_rad


def execute(block, config):

    rad, n_rad = config 
     

    #read in the model values at this sample point.
    # We look for CUDA algorithm results first; if that fails, we look for CPU-based
    # algorithm results. If that fails the module fails. This implementation assumes
    # that we are either using all CUDA algorithms or all CPU algorithms. If that is
    # not sufficient, we need a slightly more complex bit of logic here.
    try:
        #vals in datablock are saved in the following way: they are 2 dimensional, the first index is
        #the volume index, while the second one indicates the grid point. Grid is made like that:
        #it is 1-dimensional but has inside the information of the redshift z and the radii r.
        #The fastest index is r: at a fixed volume (so in a row) the first element has z_0 and r_0;
        #the second one has z_0 r_1; the third one has z_0 and r_2 ecc untill we have passed all the r index
        #Then the following element will be z_1 r_0; the next z_1 r_1 and so on until all the z values 
        #are considered. 
        profiles_miscent_model = block['SigmaMiscentY1MortCUDAScalarIntegrand', "vals"]
        profiles_cent_model = block['SigmaCentY1MortCUDAScalarIntegrand', "vals"]
        NC_miscent_model = block['NCMiscentY1MortCUDAScalarIntegrand', "vals"]
        NC_cent_model = block['NCCentY1MortCUDAScalarIntegrand', "vals"]
    except:
        profiles_miscent_model = block['SigmaMiscentY1MortScalarIntegrand', "vals"]
        profiles_cent_model = block['SigmaCentY1MortScalarIntegrand', "vals"]
        NC_miscent_model = block['NCMiscentY1MortScalarIntegrand', "vals"]
        NC_cent_model = block['NCCentY1MortScalarIntegrand', "vals"]

    miscent_fraction = block['cluster_abundance', 'fcen']

    #define the volume length from the data and the redshift one too (second index of data/length of radii)
    n_vol = profiles_miscent_model.shape[0]
    n_z = profiles_miscent_model.shape[1]/n_rad
    #change the shape of data from datablock suche that the first index is the volume, the second one
    #the redshift and the third one the radii
    profiles_miscent_model.shape = (n_vol,n_z,n_rad)
    profiles_cent_model.shape = (n_vol,n_z,n_rad)
    NC_miscent_model.shape = (n_vol,n_z,n_rad)
    NC_cent_model.shape = (n_vol,n_z,n_rad)
    miscent_fraction.shape = (n_vol,n_z,n_rad)



    '''
    #this files could be used in order to debug

    profiles_miscent_model = np.genfromtxt('sigma_y1_output/sigmamiscenty1mortscalarintegrand/vals.txt')
    profiles_cent_model = np.genfromtxt('sigma_y1_output/sigmacenty1mortScalarintegrand/vals.txt')
    NC_miscent_model = np.genfromtxt('sigma_y1_output/ncmiscentY1mortscalarintegrand/vals.txt')
    NC_cent_model = np.genfromtxt('sigma_y1_output/nccenty1mortscalarintegrand/vals.txt')
    miscent_fraction = 0.7 #np.genfromtxt('sigma_y1_output/cluster_abundance/fcen.txt') 
    '''

    #compute sigma as a 3-dimensional array with the index in the following order: volume, redshift, radii
    sigma = np.zeros(n_vol,n_z,n_rad)
    sigma = profiles_miscent_model/NC_miscent_model * miscent_fraction + (profiles_cent_model/NC_cent_model)*(1-miscent_fraction)
   


    #compute ds=deltasigma
    ds = convert_s2ds(rad, sigma)
    
    #now I change the shape of ds again in order to save it in the datablock: the first index would be the volume, while
    #the second one would be redshift*radii, with radii the faster one
    ds.shape = (n_vol, (n_z) * (n_rad))
    #now put deltasigma into the datablock.
    block["deltasigma", "deltasigma_tot"] = ds
    
    return 0
  
def cleanup(config):
    pass

  
#convert sigma in deltasigma: it returns deltasigma(called ds); profiles would be profiles=profiles(n_vol,n_z,n_rad)
def convert_s2ds(rad, profiles):

    ds = np.zeros(profiles.shape[0], profiles.shape[1], profiles.shape[2])
    #jj is the volume index,kk the redhift one and ii the radii 
    for jj in range(profiles.shape[0]):
        for kk in range(profiles.shape[1]):
            for ii in range(profiles.shape[2]):
                #compute deltasigma: np.trapz()*2/R^2 gives <sigma(<R)>, then we subtract sigma(R)
                ds_int = np.trapz(profiles[jj,kk,:(ii+1)]*rad[:(ii+1)])*2.0/(rad[ii])**2-profiles[jj,kk,ii]
                ds[jj,kk,ii] = ds_int

    return ds
    

