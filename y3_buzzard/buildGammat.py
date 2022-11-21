import os
import numpy as np
from cosmosis.datablock import option_section, names

# redshift mean with format [z0, z1, z2, z0, z1, z2, ...]
from setup_bins import zmeans_ij

#############################################
########### Build Shear Profile #############
# Estimates the gamma profiles:
#     gamma(R) = <kappa(<R)> - kappa(R)
#
# Computes the mean kappa in a radii bin.
# <kappa(R)> = 1/R \int_0^{R} kappa(R^\prime) dR^\prime
# see buildKappa.py for information on kappa
#
# Converts gamma(R) to gamma(\theta)
# Units of theta are specified on the .ini file
#
# Returns: 
# avg_shear/
# | shear_cen.txt
# | theta.txt
# | r.tx
#
# Note on the output shear shape
# for a given row there is a vector gamma(r)
# each row is a (lambda, redshift) bin 
# rows sequences are: (lbd_0, z_0), (lbd_0, z_1), ...
# r has the same number of columns of gamma(lbd, z, r)
# theta has the same shape of gamma since depends on redshift
#############################################
# Author: Johnny Esteves
# Created: 17th Nov, 2022
#############################################

cosmo = names.cosmological_parameters

# angular conversion factor
rad2deg = 180/np.pi
conv_factor = {"radians": 1.,
               "degrees": 1.*rad2deg, 
               "arcmin" : 60.*rad2deg, 
               "arcsec" : 360.*rad2deg}

def setup(options):
    section = option_section
    # radii of xi_hm, in Mpc/h
    #Mpc/h comoving distance, distance on the sky
    Radii_min = options[option_section,"Radii_min"]
    Radii_max = options[option_section,"Radii_max"]
    Radii_bins = int(options[option_section,"Radii_bins"])

    # sep_units, default is arcmin
    sep_units = str(options[option_section,"sep_units"])

    return Radii_min, Radii_max, Radii_bins, sep_units

def execute(block, config):
    Radii_min, Radii_max, Radii_bins, sep_units  = config

    # cosmological quantities
    h0 = block[cosmo, "h0"]
    z_da = block['distances', 'z']
    da = block['distances', 'd_a'] # Mpc

    # load auxialiary vectors
    r_kappa = block["avgKappaCentBu", "r"]

    # number of radial bins
    Nrbins = r_kappa.size

    # the kappa shape from the datablock is 
    # (number lambda bins, number redshift bins times number radial bins)
    # the radial bins were stride Nzbins
    kappa_cen  = block["avgKappaCentBu", "vals"]
    Nlbins, Nzr = kappa_cen.shape

    # get the number of redshift bins
    Nzbins = int(Nzr/Nrbins)

    # reshape: kappa is one profile (i.e. radial bins) for each row
    # each row is one (lambda, redshift) bin
    try:
        kappa_cen_reshaped = kappa_cen.reshape(Nlbins*Nzbins, Nrbins)
    except:
        print("Error: kappa shape was not output correctly.")
        print("Shape should (Nlbdbins, Nzbins*Nrbins")

    # interpolate on the new bining scheme
    Radii = np.logspace(np.log10(Radii_min), np.log10(Radii_max), Radii_bins)

    # compute shear profile
    # \gamma(r) = <\kappa(<r)> - k(r)
    shear_cen = np.zeros((Nlbins*Nzbin, Nrbins))
    for ij in range(Nlbins*Nzbin):
        shear_ij = compute_mean_profile(r_kappa, kappa_cen_reshaped[jj])
        shear_cen[ij] = interp1d(r, shear_ij)(Radii)
    
    # convert R [Mpc/h] to theta [arcmin/h]
    # theta depends on redshift, because of angular distance
    # for simplicity thetha will have the same shape of shear
    theta = np.zeros((Nlbins*Nzbin, Nrbins))
    for z in zmeans_ij:
        # convert R to theta
        theta[ij] = r_to_theta(Radii, z, da, z_da,
                               sep_units)

    # put into the datablock
    block["avgShear", "r"] = Radii/h0 # Mpc/h
    block["avgShear", "theta"] = theta/h0 #sep_units/h
    block["avgShear", "shear_cen"] = shear_cen
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

def r_to_theta(r, z, z_interp, da_interp, sep_units="arcmin"):
    """Physical Units to Angle
    
    Computes the angle separation in the sky
    Given a vector of angle distances and its redshift
    Interpolates and apply to the new redshift

    Note: this code is meant to recieve distances from CAMB.

    Args:
        r (array or float): physical distance separation in units of Mpc/h
        z (array or float): redshift of the object (should be same size of r)
        z_interp (array): z interpolator vector, same size of da
        da_interp (array): ang distance interpolator vector, same size of z
        sep_units (str, optional): angle units, options are `degrees`, 
        `arcmin`, `arcesec` and `radians`. Defaults to `arcmin`.

    Returns:
        theta (array): angular separation in units of sep_units
    """
    # theta = l * ang. distance [radians]
    da_sep = np.interp(z, z_da, da)
    theta_rad = r*da_sep 

    # convert to the sep_units
    theta = theta_rad*conv_factor[sep_unit]

    return theta

def cleanup(config):
    pass
