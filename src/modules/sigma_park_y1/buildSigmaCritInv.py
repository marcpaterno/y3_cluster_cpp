import os
from cosmosis.datablock import option_section, names
import numpy as np
from scipy.interpolate import interp1d
import astropy.io.fits as pyfits

#####################################################################
################### Build Sigma Critical Inverse ####################
# 
# Computes sigma critical as a function of the lens redshift
# Uses a lookup table of the effective source redshifts.
#
# In the weak lensing regime the sigma_critical is a statistical 
# measurement of the ensemble of sources. 
# 
# The effective sigma_crital can be defined as:
#     <sigma_critical(zL, zs_eff)> = 4piG/c^2 da_l <beta(zl, zs_eff)>
# Where beta is the ratio beta = dls/ds. 
# For a redshift probability distribution the mean beta
# is the average of the redshift sources.
#     <beta(zs_eff)> = \int p(zs) d(zl,zs)/d(zs) dzs
# So the effective source redshift is:
#     d(zl,zs_eff)/d(zs_eff) = <beta(zl, zs_eff)>
# which is a single number for a given lens redshift. 
###################################################################
# Author: Johnny Esteves
# Created: 17th Nov, 2022
###################################################################

cosmo = names.cosmological_parameters

def setup(options):
        section = option_section
        z_file_name=os.path.expandvars('${Y3_CLUSTER_CPP_DIR}/data/test_cluster_Y1.fits')
        z_data=pyfits.open(z_file_name)[6].data
        z_array = z_data['z_mid']

        ## FIXME DUMMY VARIABLE
        ## NEED TO CREATE A FILE WITH ZSRC_EFF
        zsrc_eff = z_data['bin1']
        
        # Load redshift bins
        z_min = options[section,"z_min"]
        z_max = options[section,"z_max"]
        z_bins = int(options[section,"z_bins"])
        zl_array = np.linspace(z_min, z_max, z_bins)

        return zsrc_eff, zl_array

def execute(block, config):
        zsrc_eff, zl_array = config

        # setup constants
        G=4.51710305e-48 #Mpc^3 / M_sol / s^2
        c=9.71561189e-15 #Mpc/s
        
        # cosmological quantities
        h0 = block[cosmo, "h0"]
        z_da = block['distances', 'z']
        da = block['distances', 'd_a'] # Mpc
        
        # interp and apply for zl, zs
        da_func = interp1d(z_da, da)
        da_zl = da_func(zl_array)
        da_zs = da_func(zsrc_eff)

        # compute sigma_crit_inv 
        sigma_crit_inv = np.zeros(len(zl_array))
        for jj, zt in enumerate(zl_array):
            da_zl_zs = da_zs - (1.0+zt)/(1.0+zsrc_eff) *da_zl[jj]
            beta = da_zl[jj] * da_zl_zs / da_zs
            sigma_crit_inv[jj] = 4.0*np.pi*G/c/c * beta
            
        block["average_sigma_crit_inv", "zlens"] = zl_array
        block["average_sigma_crit_inv", "sigma_crit"] = sigma_crit_inv
        return 0

def cleanup(config):
	pass
