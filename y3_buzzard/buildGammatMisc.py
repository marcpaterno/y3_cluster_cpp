import os
import numpy as np
from scipy.interpolate import interp1d
from cosmosis.datablock import option_section, names
from scipy import integrate

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
# | shear.txt
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
rad2deg = 180 / np.pi
conv_factor = {
    "radians": 1.,
    "degrees": 1. * rad2deg,
    "arcmin": 60. * rad2deg,
    "arcsec": 360. * rad2deg
}


def setup(options):
    section = option_section
    # radii of xi_hm, in Mpc/h
    #Mpc/h comoving distance, distance on the sky
    Radii_min = options[section, "Radii_min"]
    Radii_max = options[section, "Radii_max"]
    Radii_bins = int(options[section, "Radii_bins"])

    # sep_units, default is arcmin
    sep_units = str(options[option_section, "sep_units"])

    # grab radius bins
    r_cent = options["gammacent", "radius"]
    r_misc = options["gammamisc", "radius"]

    # do the integral over kappa
    do_int = int(options[section, "do_int"])
    return r_cent, sep_units, do_int


def execute(block, config):
    r_cent, sep_units, do_int = config

    # fraction of centered clusters
    fcen = block["cluster_abundance", "fcen"]

    # cosmological quantities
    h0 = block[cosmo, "h0"]
    z_da = block['distances', 'z']
    da = block['distances', 'd_a']  # Mpc

    # load auxialiary vectors
    # r_kappa = block["avgKappaCentBu", "radius"]

    # build average: <x> = N[x]/N[1]
    NC = block["numberCounts", "vals"]

    # number of radial bins
    Nrbins = r_cent.size

    # the kappa shape from the datablock is
    # (number lambda bins, number redshift bins times number radial bins)
    # the radial bins were stride Nzbins
    gamma_cen = block["gammacent", "vals"]
    gamma_mis = block["gammamisc", "vals"]
    Nlbins, Nzr = gamma_cen.shape

    # get the number of redshift bins
    Nzbins = int(Nzr / Nrbins)
    # print("Number of redshift bins: %i" % Nzbins)

    # reshape: kappa is one profile (i.e. radial bins) for each row
    # each row is one (lambda, redshift) bin
    try:
        gamma_cen_reshaped = gamma_cen.reshape(Nlbins * Nzbins, Nrbins)
        gamma_mis_reshaped = gamma_mis.reshape(Nlbins * Nzbins, Nrbins)
    except:
        print("Error: kappa shape was not output correctly.")
        print("Shape should (Nlbdbins, Nzbins*Nrbins")

    # compute shear profile
    # \gamma(r) = <\kappa(<r)> - k(r)
    shear_cen = np.zeros((Nlbins * Nzbins, Nrbins))
    shear_mis = np.zeros((Nlbins * Nzbins, Nrbins))
    for ij in range(Nlbins * Nzbins):
        shear_cen[ij] = gamma_cen_reshaped[ij]/NC.flatten()[ij]
        shear_mis[ij] = gamma_mis_reshaped[ij]/NC.flatten()[ij]

    # convert R [Mpc/h] to theta [arcmin/h]
    # theta depends on redshift, because of angular distance
    # for simplicity theta will have the same shape of shear
    theta = np.zeros((Nlbins * Nzbins, Nrbins))
    for ij, z in enumerate(zmeans_ij):
        # convert R to theta
        theta[ij] = r_to_theta(r_cent, z, z_da, da, sep_units)

    # put into the datablock
    block["shear", "r"] = r_cent / h0  # Mpc/h
    block["shear", "theta"] = theta / h0  #sep_units/h
    block["shear", "shear_cent"] = shear_cen
    block["shear", "shear_misc"] = shear_mis
    block["shear", "shear_full"] = fcen*shear_cen + (1-fcen)*shear_mis
    return 0


def compute_mean_profile(r, fx, rmin=0.1):
    """Computes \Delta f(x) = <f(<x)> - f(x)
    
    Average profiles excess of f(x) 
    Example: f(x) = \Sigma
    """
    # start the integration
    # delta_profile = np.zeros(r.size)
    # for ii, ri in enumerate(r[1:]):
    #     ix = np.where(r<ri)[0]
    #     #delta_profile[ii]  = np.trapz(fx[ix], x=r[ix])/ri - fx[ii]
    profile = interp1d(r, fx, fill_value='extrapolate')

    delta_profile = np.full(r.size, np.nan)
    for ii, ri in enumerate(r):
        if ri > rmin:
            delta_profile[ii] = integrate.quad(profile, 0.05,
                                               ri)[0] / ri - profile(ri)
            #delta_profile[ii] = integrate.fixed_quad(profile, 0.001, ri)[0]/ri - profile(ri)
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
    da_sep = np.interp(z, z_interp, da_interp)
    theta_rad = (r / da_sep)

    # convert to the sep_units
    theta = theta_rad * conv_factor[sep_units]

    return theta


def cleanup(config):
    pass
