import os
import numpy as np
from scipy.interpolate import interp1d
from cosmosis.datablock import option_section, names

# FIDUCIAL COSMOLOGY
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)

#####################################################################
###################### Scale Shift Cosmology  #######################
# There are some distances that we need to define a fiducial cosmology
# To adapt to a given cosmology we can re-scale the distances.
#
# scaleShiftCosmo is basically the ratio of the comoving distance
# and the hubble factor
#     scaleShiftCosmo = (dist_c/h) / (dist_c_fid/h_fid)
#
# our fiducial cosmology is set as Omega_m = 0.3 and H0 = 70
###################################################################
# Author: Johnny Esteves
# Created: 18th Nov, 2022
###################################################################

cosmo = names.cosmological_parameters

def setup(options):
        section = option_section

        return 0

def execute(block, config):        
        # cosmological quantities
        h0 = block[cosmo, "h0"]
        z_dc = block["distances",'z']
        dc = block["distances",'d_a']*(1+z_dc) # comoving distance Mpc

        # fiducical cosmology
        h0_fid = cosmo.H0/100.
        dc_fid = cosmo.comoving_distance(z_dc).value

        # scale shift
        scale_shift = (dc/h0)/(dc_fid/h0_fid)

        block["scaleShiftCosmo", "z"] = z_dc
        block["scaleShiftCosmo", "sigma_crit"] = scale_shift
        return 0

def cleanup(config):
	pass
