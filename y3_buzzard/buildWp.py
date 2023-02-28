# -*- coding: utf-8 -*-
import os
import numpy as np
from scipy.interpolate import interp1d
from cosmosis.datablock import option_section, names
from scipy import integrate

# redshift mean with format [z0, z1, z2, z0, z1, z2, ...]
from setup_bins import zmeans_ij

#################################################
########### Build Wp(\thea) Profile #############
# Projects \xi(r) onto the sky plane:
#     Wp(R) = \int \xi(\sqrt{R^2+\pi^2}) d\pi
#
# Which is the Abell tranformation.
# The limits of the integration are -\pi_max, \pi_max
#
# Converts Wp(R) to Wp(\theta) and corrects for the
# fiducial cosmology.
#
# Returns:
# wp/
# | wp_cen.txt
# | theta.txt
# | r.tx
#
# Note on the output wp shape
# for a given row there is a vector wp(r)
# each row is a (lambda, redshift) bin
# rows sequences are: (lbd_0, z_0), (lbd_0, z_1), ...
# r has the same number of columns of wp(lbd, z, r)
# theta has the same shape of wp since depends on redshift
#############################################
# Author: Johnny Esteves
# Created: 10th Dec, 2022
#
# Based on buildGammat.py
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

    # pimax; radial parallel component max value
    pimax = float(options[option_section, "pimax"])

    # grab radius bins
    r_xi = options["avgWp", "radius"]

    # do the integral over kappa
    do_int = int(options[section, "do_int"])

    return Radii_min, Radii_max, Radii_bins, r_xi, sep_units, do_int


def execute(block, config):
    Radii_min, Radii_max, Radii_bins, r_xi, sep_units, do_int = config

    # cosmological quantities
    h0 = block[cosmo, "h0"]
    z_da = block['distances', 'z']
    da = block['distances', 'd_a']  # Mpc

    # build average: <x> = N[x]/N[1]
    NC = block["numberCounts", "vals"]

    # number of radial bins
    Nrbins = r_xi.size

    # the \xi shape from the datablock is
    # (number lambda bins, number redshift bins times number radial bins)
    # the radial bins were stride Nzbins
    wp_cen = block["avgWp", "vals"]
    Nlbins, Nzr = wp_cen.shape

    # get the number of redshift bins
    Nzbins = int(Nzr / Nrbins)
    print("Number of redshift bins: %i" % Nzbins)

    # reshape: wp is one profile (i.e. radial bins) for each row
    # each row is one (lambda, redshift) bin
    try:
        wp_cen_reshaped = wp_cen.reshape(Nlbins * Nzbins, Nrbins)
    except:
        print("Error: kappa shape was not output correctly.")
        print("Shape should (Nlbdbins, Nzbins*Nrbins")

    # interpolate on the new bining scheme
    Radii = np.logspace(np.log10(Radii_min), np.log10(Radii_max), Radii_bins)

    # compute the abell transform
    # Wp(R) = \int \xi(\sqrt{R^2+\pi^2}) d\pi
    wp_out = np.zeros((Nlbins * Nzbins, Radii_bins))
    for ij in range(Nlbins * Nzbins):
        if do_int:
            wp_ij = compute_abell_transform(r_xi, wp_cen_reshaped[ij], pimax)
        else:
            wp_ij = wp_cen_reshaped[ij]

        wp_avg_ij = wp_ij / NC.flatten()[ij]
        wp_out[ij] = interp1d(r_xi, wp_avg_ij, bounds_error=False)(Radii)

    # convert R [Mpc/h] to theta [arcmin/h]
    # theta depends on redshift, because of angular distance
    # for simplicity thetha will have the same shape of shear
    theta = np.zeros((Nlbins * Nzbins, Radii_bins))
    for ij, z in enumerate(zmeans_ij):
        # convert R to theta
        theta[ij] = r_to_theta(Radii, z, z_da, da, sep_units)

    # put into the datablock
    block["wp", "r"] = Radii / h0  # Mpc/h
    block["wp", "theta"] = theta / h0  #sep_units/h
    block["wp", "wp_cen"] = wp_out
    return 0


def compute_abell_transform(r, fr, pimax=100., int_func=np.trapz):
    # Computes Wp(R) = \int \xi(\sqrt{R^2+\pi^2}) d\pi
    #
    # The equation above can be re-written to:
    #
    # Wp(R) = \int (dr/\sqrt{r^2-R^2}) r \xi(r).
    #
    # Calculation of the integral  used in Abel transform
    #          r_max
    #         ⌠
    #         ⎮     r \xi(r)
    #         ⎮ ────────────── dr
    #         ⎮    ___________
    #         ⎮   ╱  2   2
    #         ⎮ ╲╱  r - R
    #         ⌡
    #         R
    # Where r_max=\sqrt{R^2+\pi_{max}^2}
    # Code based on: PyAbel package, direct method.
    # https://github.com/PyAbel/PyAbel/blob/master/abel/direct.py
    # TODO: Implement the C++ version.

    # Switch to the r coordinates.
    fr *= 2 * r

    if is_uniform_sampling(r):
        int_opts = {'dx': abs(r[1] - r[0])}
    else:
        int_opts = {'x': r}

    # Initializiate output
    out = np.zeros(fr.size)

    # Initialize grid
    R, rr = np.meshgrid(r, r, indexing='ij')
    i_vect = np.arange(len(r), dtype=int)
    II, JJ = np.meshgrid(i_vect, i_vect, indexing='ij')

    # Set integration limits
    # only for R < r <= rmax
    mask = (II < JJ)
    mask &= (rr <= np.sqrt(R**2 + pimax**2))

    I_sqrt = np.zeros(R.shape)
    I_sqrt[mask] = np.sqrt((rr**2 - R**2)[mask])

    I_isqrt = np.zeros(R.shape)
    I_isqrt[mask] = 1. / I_sqrt[mask]

    # create a mask that just shows the first two points of the integral
    mask2 = ((II > JJ - 2) & (II < JJ + 1))

    # perfom the integration
    Integrand = fr[None, :] * I_isqrt  # set up the integral
    out = int_func(Integrand, axis=1, **int_opts)  # take the integral
    # correct for the extra triangle at the start of the integral
    out = out - 0.5 * int_func(Integrand * mask2, axis=1, **int_opts)

    return out


def is_uniform_sampling(r):
    """
    Returns True if the array is uniformly spaced to within 1e-13.
    Otherwise False.
    """
    dr = np.diff(r)
    ddr = np.diff(dr)
    return np.allclose(ddr, 0, atol=1e-13)


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
