import os
import numpy as np
from cosmosis.datablock import option_section, names

from scipy.interpolate import interp1d
from scipy.special import erf
import cluster_toolkit as ct

from astropy.constants import G

from haloModel import biasModel, lensingModel, scaleShiftCosmo

##################################################################
############# Build Wp(R,z) and Gammat(R,z) ############
# Estimates Halo Model Parameters
# Wp, \Sigma, \Delta \Sigma and bias(z)
# 
#############################################
# Author: Johnny Esteves
# Created: May 16, 2024
#############################################

import astropy.cosmology
cosmo_names = names.cosmological_parameters

def setup(options):
    section = option_section
    #Mpc/h comoving distance, distance on the sky
    R_perp_min = options[section,"R_perp_min"]
    R_perp_max = options[section,"R_perp_max"]
    R_perp_bins = int(options[section,"R_perp_bins"])

    # radii of xi_hm, in Mpc/h
    Radii_min = options[section,"Radii_min"]
    Radii_max = options[section,"Radii_max"]
    Radii_bins = int(options[section,"Radii_bins"])

    #mass (float): Halo mass Msun/h.
    M_min = options[section,"M_min"]
    M_max = options[section,"M_max"]
    M_bins = int(options[section,"M_bins"])

    params_out = R_perp_min, R_perp_max, R_perp_bins, \
    Radii_min, Radii_max, Radii_bins, M_min, M_max, M_bins
    return params_out
    

def execute(block, config):
    section_name = "haloModel"
    
    R_perp_min, R_perp_max, R_perp_bins, Radii_min, Radii_max,\
    Radii_bins, M_min, M_max, M_bins  = config

    # cosmo parameters
    omega_m = block[cosmo_names, "omega_m"]
    omega_b = block[cosmo_names, "omega_b"]
    H0 = block[cosmo_names, "H0"]

    cosmology = astropy.cosmology.FlatLambdaCDM(H0*100, omega_m, Ob0=omega_b, Tcmb0=2.725)

    # load camb matter-matter power spectrums
    # linear is used for the bias
    k_h = block["matter_power_lin", "k_h"]
    P_k = block["matter_power_lin", "p_k"]
    z_k = block["matter_power_lin", "z"]

    P_k_nl = block["matter_power_nl", "p_k"]
    k_nl = block["matter_power_nl", "k_h"]
    z_nl = block["matter_power_nl", "z"]

    z = z_k

    z_az = block["growth_parameters", "z"]
    daz = block["growth_parameters", "d_z"]
    daz = np.interp(z, z_az, daz)

    # compute overdensities; physical not comoving
    rho_c0 = float(cosmology.critical_density0.to('Msun/Mpc^3').value)
    # rho_cz = cosmology.critical_density(z).to('Msun/Mpc^3').value
    rho_m = omega_m*rho_c0
    rho_mz = rho_m*(1.+z)**3
    
    # setup bins
    nz = len(z)
    R_perp = np.logspace(np.log10(R_perp_min), np.log10(R_perp_max), R_perp_bins)
    Radii = np.logspace(np.log10(Radii_min), np.log10(Radii_max), Radii_bins)
    M = np.logspace(np.log10(M_min), np.log10(M_max), M_bins)
    logM = np.log(M)

    # Compute Bias (M, Z)
    # The bias is computed at the peak height (z=0) from cluster toolkit
    # The peak-height evolution is computed from the growth function D(z)
    # The vector Bias has shape (z.size, M.size)
    bM = biasModel(k_h, P_k[0], omega_m, odelta=200)
    nu = bM.nu_at_M(M)
    Bias = np.array([bM.bias_at_nu(nu/dazi) for dazi in daz])

    #### Step 2) Sigma, deltaSigma
    # Initialize the lensing module
    lensModel = lensingModel(R_perp, omega_m, 200)

    # First halo term is a NFW profile with mass M and concentration based on relation Child18
    # Second halo term is the halo-halo correlation function projected along the line of sight (i.e. missing the halo bias)
    lensModel.first_halo_term(M, z=0, conc_model_name="Child18")
    lensModel.second_halo_term(z, k_nl, P_k_nl)

    # Step 5) shift cosmological scale
    # used to shift gammat to the fiducial cosmology used in the datavecotr measurement
    scaleShift, hubbleShift = scaleShiftCosmo(z, cosmology)

    # put into the datablock
    block[section_name, "m_h"] = M
    block[section_name, "lnM"] = logM
    block[section_name, "z"] = z
    block[section_name, "rhoc"] = rho_mz

    # Wp
    block[section_name, "Rp"] = Radii
    block[section_name, "Wp_hh"] = lensModel.Wp

    # Sigma, deltaSigma [Msun/pc^2]
    block[section_name, "r_sigma"] = R_perp
    block[section_name, "Sigma_hh"] = lensModel.Sigma['2h'] 
    block[section_name, "Sigma_nfw"] = lensModel.Sigma['1h'] 
    block[section_name, "dSigma_hh"] = lensModel.dSigma['2h'] 
    block[section_name, "dSigma_nfw"] = lensModel.dSigma['1h'] 
    block[section_name, "concentration"] = lensModel.c
    
    # Bias
    block[section_name, "bias"] = Bias

    # damped Pk
    block[section_name, "scale_shift"] = scaleShift
    block[section_name, "hubble_shift"] = hubbleShift
    block[section_name, "k"] = k_h

    return 0

def cleanup(config):
    pass