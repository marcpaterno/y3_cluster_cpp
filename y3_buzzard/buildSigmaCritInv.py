import os
import numpy as np
from scipy.interpolate import interp1d
from cosmosis.datablock import option_section, names

# redshift mean with format [z0, z1, z2]
from setup_bins import zmeans

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
#     d(zl,zs_eff)/d(zs_eff) = <beta(zl, zs)>
# which is a single number for a given lens redshift. 
#
#####################################################################
#################### Code Implementation ############################
# We import beta from a lookup table. 
# Since the quantity is cosmological we need to apply a cosmological
# shift on the distance. This correction is:
#     beta = 1 - scale_ratio(zs, zl) * (1-beta_fid)
#
# where beta_fid is the one computed by the fiducial cosmology and
#     scale_ratio(zs, zl) = scaleShiftCosmo(zs)/scaleShiftCosmo(zl)
#
# the scaleShiftCosmo is basically the ratio of the comoving distance
# and the hubble factor
#     scaleShiftCosmo = (dist_c/h) / (dist_c_fid/h_fid)
#
# for each cosmological step we compute the cosmological shift 
# on the scaleShiftCosmo.py module
###################################################################
# Author: Johnny Esteves
# Created: 17th Nov, 2022
###################################################################

cosmo = names.cosmological_parameters

def setup(options):
        section = option_section

        # loads beta table
        beta_lut_fname = os.path.expandvars('${Y3_CLUSTER_CPP_DIR}/data/beta_buzzard_interp.npz')
        betaTable = np.load(beta_lut_fname)

        # interpolates beta over lens redshift
        beta_fid = interp1d(betaTable['z'], betaTable['beta'])(zl_array)
        z_src = betaTable['zsrc']

        return z_src, beta

def execute(block, config):
        z_src, beta_fid_zl  = config

        # lenses redshift is equal to the cluster mean redshifts
        z_lenses = zmeans

        # setup constants
        G=4.51710305e-48 # Mpc^3 / M_sol / s^2
        c=9.71561189e-15 # Mpc/s
        
        # cosmological quantities
        h0 = block[cosmo, "h0"]
        z_da = block["distances",'z']
        da = block["distances",'da'] # Mpc

        # lens comoving distance
        dc_lenses = interp1d(z_da, da)(z_lenses) * (1+z_lenses)

        # cosmological shift := (H/Hfid)*(dc/dc_fid)
        # fid stands for the fiducial cosmology
        z_shift = block['cosmoShift', 'z']
        shift = block['cosmoShift', 'shift']

        # compute the shift ratio of dls/ds
        scale_shift_func = interp1d(z_shift, shift)
        scale_ratio = scale_shift_func(z_src)/scale_shift_func(z_lenses)

        # beta with the cosmological correction
        # beta = dls/ds is unitless
        beta = 1 - scale_ratio*(1-beta_fid)
        # if scale_ratio is one, beta equals to beta_fid

        # compute sigma_crit_inv [Msun]
        sigma_crit_inv = 4.0*np.pi*G/c/c * dc_l * beta
        
        block["sigma_crit_inv", "z"] = z_lenses
        block["sigma_crit_inv", "sigma_crit"] = sigma_crit_inv
        return 0

def cleanup(config):
	pass
