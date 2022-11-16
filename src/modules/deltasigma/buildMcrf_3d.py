import os
from cosmosis.datablock import option_section, names
import numpy as np
from scipy.special import erf
import cluster_toolkit as ct

# buildMcrf_3d = buld matter-matter- 3d-correlation function
# but damped as in Estrada et al 2009  eq 8
# https://ui.adsabs.harvard.edu/abs/2009ApJ...692..265E/abstract

cosmo = names.cosmological_parameters


def setup(options):

    # radii of xi_hm, in Mpc/h
    #Mpc/h comoving distance, distance on the sky
    Radii_min = options[option_section,"Radii_min"]
    Radii_max = options[option_section,"Radii_max"]
    Radii_bins = int(options[option_section,"Radii_bins"])
    sigma_photoz= options[option_section,"sigma_photoz"]

    return Radii_min, Radii_max, Radii_bins, sigma_photoz

def execute(block, config):
    Radii_min, Radii_max, Radii_bins, sigma_photoz  = config

    k = block["matter_power_lin", "k_h"]
    P_k = block["matter_power_lin", "p_k"]
    z_k = block["matter_power_lin", "z"]

    k_nl = block["matter_power_nl", "k_h"]
    P_k_nl = block["matter_power_nl", "p_k"]
    z_k_nl = block["matter_power_nl", "z"]


    Radii = np.logspace(np.log10(Radii_min), np.log10(Radii_max), Radii_bins)
    dampingFactor = np.pi*erf(sigma_photoz * k_nl)/(2*sigma_photoz*k_nl)

    nz = z_k_nl
    Xi_mm = np.zeros((nz, Radii_bins))
    for i in range(nz):
        damped_Pk_nl[i]  = dampingFactor * P_k_nl[i]
        xi_mm_at_z = ct.xi.xi_mm_at_r(Radii, k_nl, damped_P_k_nl[i])
        Xi_mm[i] = xi_mm_at_z

    block["correlationFunction", "Damped_Xi_mm"] = Xi_mm
    block["correlationFunction", "R_Xi"] = Radii
    block["correlationFunction", "z"] = z
    return 0


def cleanup(config):
    pass
