from __future__ import print_function
from cosmosis.datablock import names
from cosmosis.datablock import option_section
import numpy as np
import scipy.integrate
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d

cosmo = names.cosmological_parameters

# Jim Annis Nov 2022
#
# Code is descended from 
#   cosmosis-standard-library/bolztmann/sigmar/sigmar.py
# modified to calculate the collapse mass scale following
# Child, Habib, Heitmann et al 2018, eq 13

# Option setup part.  Read options from in ifile.
# Definition of z-bins and R-bins
# Specify which camb power spectrum: 
#   matter_power = matter_power_lin, usually


def setup(options):
    zmin = options[option_section, "zmin"]
    zmax = options[option_section, "zmax"]
    dz = options[option_section, "dz"]
    z = np.arange(zmin, zmax, dz)
    rmin = options[option_section, "rmin"]
    rmax = options[option_section, "rmax"]
    dr = options[option_section, "dr"]
    R = np.arange(rmin, rmax, dr)
    crop_klim = options.get_bool(option_section, "crop_klim", True)
    R = np.atleast_1d(R)
    z = np.atleast_1d(z)
    blockname = options[option_section, "matter_power"]
#    print("Sigmar(R,z) will be evaluated at:")
#    print("z = ", z)
#    print("R = ", R)
#    print("crop_klim = ", crop_klim)
    return (z, R, blockname, crop_klim)


def sigint(lnk, r, z, rbs):
    k = np.exp(lnk)
    x = k * r
    w = 3 * (-x * np.cos(x) + np.sin(x)) / x**3
    p = rbs.ev(k, z)
    tmp = w**2 * k**3 * p / (2 * np.pi**2)
    return tmp


def execute(block, config):
    z, R, blockname, crop_klim = config

    karray, zarray, powerarray = block.get_grid(blockname, "k_h", "z", "p_k")
    omega_m = block[cosmo, "omega_m"]
    h0 = block[cosmo, "h0"]

    rbs = RectBivariateSpline(karray, zarray, powerarray)
    kmin_overall = np.log(karray.min())
    kmax_overall = np.log(karray.max())

    kmin = kmin_overall
    kmax = kmax_overall

    sigma2r = 0

    sigma2r = np.zeros((np.size(R), np.size(z)))
    for i, rloop in enumerate(R):

        if crop_klim:

            kmin = max(np.log(.01 / rloop), kmin_overall)
            kmax = min(np.log(100. / rloop), kmax_overall)

        for j, zloop in enumerate(z):

            sigma2r[i, j] = scipy.integrate.quad(
                sigint, kmin, kmax, args=(rloop, zloop, rbs), epsrel=1e-6)[0]


    Rstar = np.ones(z.size)
    interp2 = interp2d(z,R,sigma2r)
    for i in range(z.size):
        zed = z[i]
        s2 = interp2(zed,R) # returns a vector of sigma(r) at z
        s2 = np.squeeze(s2)
        # Johnny on 20th Nov added fill_value="extrapolate"
        interp1 = interp1d(s2, R, fill_value="extrapolate")  # returns 1d interp object returning r given sigma
        Rstar[i] = interp1(1.686**2)                # collapse radial scale, Child's eq 13

    rhocrit = 2.77533742639e+11  #h^2 M_sun Mpc^-3
    rhocrit_z = rhocrit*(1+z)**3
    Mstar = (4.*np.pi/3.) *rhocrit * omega_m*h0**2  * Rstar**3

    block.put("mstar", "z", z)
    block.put("mstar", "mstar", Mstar)
    return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness
    return 0
