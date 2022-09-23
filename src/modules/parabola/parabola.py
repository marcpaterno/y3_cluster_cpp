import os
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section, names
import numpy as np
from scipy.interpolate import interp1d, interp2d
from scipy.interpolate import RectBivariateSpline
import scipy.special
import cluster_toolkit as ct
from scipy import optimize
import matplotlib.pyplot as plt



def setup(options):
        
    section = option_section
    
    #radii
    rad = options[section, "radius"]

    #length of radii
    n_rad = len(rad)
    #redshift and richness values
    z_mean = option[section, "z_mean"]
    l_mean = option[section, "l_mean"]

    return rad, n_rad, z_mean, l_mean


def execute(block, config):
    rad, n_rad, z_mean, l_mean = config
        

    #load deltasigma_tot from finish.py
    #now put deltasigma into the datablock. It is 2 dimensional, with the first index which is len(volume)
    #and the second one is len(rad)*len(redshift), with radii_index the faster one
    ds =  block["deltasigma", "deltasigma_tot"]

    #from ds we can understand len(volume) and len(z)
    n_vol = ds.shape[0]
    n_z = ds.shape[1]/n_rad

    #I put a different shape for ds, so that the first index is volume,the second redshift and the third radii
    ds.shape = (n_vol,n_z,n_rad)

    #load the parabola parameters, at each step of the MCMC they are constants: they are shared with all bins
    a0 = block["parabola", "a_0"] 
    alpha_a = block["parabola", "alpha_a"] 
    beta_a = block["parabola", "beta_a"] 
    b0 = block["parabola", "b_0"] 
    alpha_b = block["parabola", "alpha_b"] 
    beta_b = block["parabola", "beta_b"]
    c0 = block["parabola", "c_0"] 
    alpha_c = block["parabola", "alpha_c"]  
    beta_c = block["parabola", "beta_c"]

        
        
    #model for the parabola function, theta:parameters, x will be logM, z:redshift, l:richness
    def model(theta, x, z, l):
        a0, alpha_a, beta_a, b0, alpha_b, beta_b, c0, alpha_c, beta_c = theta
        aa = a0 + ( l / 30 )**(alpha_a) + ((1+z)/1.3)**(beta_a)
        bb = b0 + ( l / 30 )**(alpha_b) + ((1+z)/1.3)**(beta_b)
        cc = c0 + ( l / 30 )**(alpha_c) + ((1+z)/1.3)**(beta_c)
        model = aa * x**2 + bb * x + cc
        return model
        
    #define theta:parabola parameters list, ds_new:the old deltasigma(ds) multiplied for the parabola
    theta = a0, alpha_a, beta_a, b0, alpha_b, beta_b, c0, alpha_c, beta_c 
    

    #define the new deltasigma, which will be the old one multiplied by parabola
    ds_new = np.zeros(n_vol,n_z,n_rad)
        
    #multiply deltasigma for the parabola
    for ii in range(n_vol):
        for jj in range(n_z):
            deltasigma = ds[ii,jj,:] * (model (theta, rad, z_mean[jj], l_mean[ii])) 
            ds_new[ii,jj,:] = deltasigma
            
        
    #overwrite in the datablock the new data, deltasigma is the old deltasigma multiplied for the parabola
    #ds_new is 2-dimensional: the first index is volume_index;
    #the second index is redshift_index*radii_index, and the faster is the radii index.
    #we need to save data in datablock with just 2 dimensional things in order to run cosmosis,
    #that is why we need to create a sort of grid for the second index
    ds_new.shape = (n_vol,(n_z)*(n_rad))
    block["deltasigma", "deltasigma_tot"] = ds_new
    return 0


def cleanup(config):
        pass
