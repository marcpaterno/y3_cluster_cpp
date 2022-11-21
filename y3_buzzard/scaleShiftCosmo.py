import os
import numpy as np
from scipy.interpolate import interp1d
from cosmosis.datablock import option_section, names

# FIDUCIAL COSMOLOGY
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)

#####################################################################
###################### Scale Shift Cosmology  #######################
# In some parts of the code, we have quantities that are cosmology
# dependent on the comoving/angular distance. 
# To adapt to a new cosmology we can re-scale the distances by taking
# into account the fiducial cosmolgoy.
#
# scaleShiftCosmo is basically the ratio of the comoving distance
# and the hubble factor
#     scaleShiftCosmo = (dist_c/h) / (dist_c_fid/h_fid)
#
# the fiducial cosmology was set to Omega_m = 0.3 and H0 = 70
###################################################################
# Author: Johnny Esteves
# Created: 18th Nov, 2022
###################################################################

def scaleShiftCosmo(znew, block, h0=0.7):
        # cosmological quantities
        # h0 = block[cosmo, "h0"]

        z_dc = block["distances",'z']
        dc = block["distances",'d_a']*(1+z_dc) # comoving distance Mpc

        # fiducical cosmology
        h0_fid = cosmo.H0/100.
        dc_fid = cosmo.comoving_distance(z_dc).value

        # scale shift
        scale_shift_vec = dc/dc_fid

        # interpolate for the new redshift
        scale_shift = np.interp(znew, z_dc, scale_shift_vec)
        return scale_shift