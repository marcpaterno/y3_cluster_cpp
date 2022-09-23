import numpy as np
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.interpolate import interp1d

##The aim of this module is to convert sigma into deltasigma
#this is made up with datablock structure, without changing the shape 
def setup(options):
   
    section = option_section
    
   
    # model module     
    rad = options[section, "radius"]
    #Redges = np.logspace(np.log10(0.0323), np.log10(30.), num=15+1)
    
    
    return rad
  
  
  
  
def execute(block, config):
  
    
    rad = config
     

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

    '''
    #this files could be used in order to debug

    profiles_miscent_model = np.genfromtxt('sigma_y1_output/sigmamiscenty1mortscalarintegrand/vals.txt')
    profiles_cent_model = np.genfromtxt('sigma_y1_output/sigmacenty1mortScalarintegrand/vals.txt')
    NC_miscent_model = np.genfromtxt('sigma_y1_output/ncmiscentY1mortscalarintegrand/vals.txt')
    NC_cent_model = np.genfromtxt('sigma_y1_output/nccenty1mortscalarintegrand/vals.txt')
    miscent_fraction = 0.7 #np.genfromtxt('sigma_y1_output/cluster_abundance/fcen.txt') 
    '''

    #compute sigma  
    sigma = profiles_miscent_model/NC_miscent_model * miscent_fraction + (profiles_cent_model/NC_cent_model)*(1-miscent_fraction)
   


    #compute ds=deltasigma and aver_ds=average deltasigma    
    #ds, aver_ds = convert_s2ds(rad, sigma, Redges)
    ds = convert_s2ds(rad, sigma)
    
    #make aver_ds into one dimension
    #model_vector = aver_ds.flatten()
    
    #now put deltasigma into the datablock. It is 2 dimensional, with the first index the volume one
    #and the second redshift_index * radii_index
    block["deltasigma", "deltasigma_tot"] = ds
    
    return 0
  
def cleanup(config):
    pass

  
#convert sigma in deltasigma: it returns deltasigma(called ds) and the average deltasigma
def convert_s2ds(rad, profiles):
    #print(ds.shape, profiles.shape)
    #profiles.shape[1] is the lenght of a row, so (z_bins)*(r_bins) so ngrid=z_bins
    ngrid=int(profiles.shape[1]/len(rad))
    #compute ds=ds(len(volume)*len(z) , len(r))
    ds= np.zeros([profiles.shape[0]*ngrid, len(rad)])
    #jj represents the volume index; kk represents the z (redshift) index; ii is the radii index
    for jj in range(profiles.shape[0]):
        for kk in range(ngrid):
            for ii in range(len(rad)):
                #profiles_ind is z*r_index with r the fastest index: if we are considering z_0 (so kk=0)
                #then profile_ind=0; if we consider z_1 (kk=1) then profiles_ind=len(r); if z_2 then
                #profile_ind=2*len(r). So profile_ind indicates the first element in the grid which has 
                #a new value of the redshift, different from the previous element. It is the beginning 
                #point for a fixed z and after that we have len(r) elements with the same z but different r
                profiles_ind = kk * len(rad)
                #profiles[jj,profiles_ind:(profiles_ind+ii+1)] are the profile values with jj-volume,
                #at a fixed redshift (the kk one) considered from the first value of r (r[0]) to 
                #the radii value r[ii] (r[0] and r[ii] are included). These values are then respectively
                #multiplied for the radii values in range between r[0] and r[ii] (included). That quantities 
                #when we will use sigma as profiles, are just Sigma(r)*r evaluated at different r. That is our 
                #y, so the input array that we are integrating.
                #the next argument of the function represents the sample points corresponding to the y values. 
                #We consider all the radii between r[0] and r[ii] (included), multiply them for 2 and then divide
                #for r[ii]**2. That is our x points; we are also considering the normalization 2/r**2, so doing the 
                #integral from 0 to r[ii] of sigma(r)*r in dr, and multiplying it for 2/r[ii]**2 
                #we obtain the average sigma that we have for radii<r[ii], so it is our <sigma(<r[ii])>.
                #Then we subtract profiles[jj, [profiles_ind+ii)] which is just sigma(r[ii]) at that
                #redshift (kk) and for that volume (jj)
                ds_int = np.trapz(profiles[jj, profiles_ind:(profiles_ind+ii+1)] *rad[:(ii+1)],
                    rad[:(ii+1)]) *2.0/(rad[ii])**2 
                    - profiles[jj, (profiles_ind+ii)]
                #save the data in the same way of the datablock
                ds[jj,profiles_ind+ ii] = ds_int
    '''
    #compute res which is the average deltasigma
    #res(len(z)*len(volume), len(Redges)-1)
    res = np.zeros([ds.shape[0], len(Redges)-1])
    #jj is the index that contains information about z and vol
    for jj in range(ds.shape[0]):
        #we interpolate ds with the radii
        prof_interp=interp1d(rad, ds[jj, :])
        #prof_interp=interp1d(rad, profiles[jj, :])
        #we define a function: it gives prof_interp(x) * x, so deltasigma(x)*x
        prof_func=lambda x: prof_interp(x) * x
        for ii in range(len(Redges)-1):
            # averaging DS ranges
            #we integrate prof_func between Redges[ii] and Redges[ii+1] (the last is from
            #Redges[len(Redges)-1] to Redges[len(Redges)]. res_int is the integral with unc_int is an estimate of
            #the absolute error in the result.
            res_int, unc_int = integrate.quad(prof_func, Redges[ii], Redges[ii+1])
            res_int = res_int *2/((Redges[ii+1]**2) - (Redges[ii]**2))
            #res is the average deltasigma that we have between Redges[ii] and Redges[ii+1] 
            #print(jj, ii, res_int)
            res[jj, ii] = res_int
    '''        
    return ds
