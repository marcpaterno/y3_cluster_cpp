import os
from cosmosis.datablock import option_section, names
import numpy as np
from scipy.special import erf
import cluster_toolkit as ct

#############################################
########### Build Shear Profile #############
# Computes the cumulative kappa in a bin
# gamma = <kappa>(<r) - kappa(r)
#############################################
# Author: Johnny Esteves
# Created: 17th Nov, 2022
#############################################

cosmo = names.cosmological_parameters

def setup(options):
    section = option_section
    # radii of xi_hm, in Mpc/h
    #Mpc/h comoving distance, distance on the sky
    Radii_min = options[option_section,"Radii_min"]
    Radii_max = options[option_section,"Radii_max"]
    Radii_bins = int(options[option_section,"Radii_bins"])

    return Radii_min, Radii_max, Radii_bins, sigma_photoz

def execute(block, config):
    Radii_min, Radii_max, Radii_bins, sigma_photoz  = config
    # load kappa
    r = block["avg_kappa", "r"]
    kappa_cen  = block["avg_kappa", "kappa_cen"]
    kappa_mis = block["avg_kappa", "kappa_mis"]
    
    # compute shear profile
    # \gamma(r) = <\kappa(<r)> - k(r)
    shear_cen = compute_mean_profile(r, kappa_cen)
    shear_mis = compute_mean_profile(r, kappa_mis)
    
    # interpolate on the new bining scheme
    Radii = np.logspace(np.log10(Radii_min), np.log10(Radii_max), Radii_bins)
    avg_shear_cen = interp1d(r, shear_cen)(Radii)
    avg_shear_mis = interp1d(r, shear_mis)(Radii)
    
    # put into the datablock
    block["avg_shear", "r"] = Radii
    block["avg_shear", "shear_cen"] = avg_shear_cen
    block["avg_shear", "shear_mis"] = avg_shear_mis
    return 0

def compute_mean_profile(r, fx):
    """Computes \Delta f(x) = <f(<x)> - f(x)
    
    Average profiles excess of f(x) 
    Example: f(x) = \Sigma
    """
    # start the integration
    delta_profile = np.zeros(r.size)
    for ii, ri in enumerate(r):
        ix = np.where(r<ri)[0]
        delta_profile[ii]  = np.trapz(fx[ix], x=r[ix])/ri - fx[i]
    return delta_profile    

def cleanup(config):
    pass
