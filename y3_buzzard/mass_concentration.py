from __future__ import print_function
from cosmosis.datablock import names
from cosmosis.datablock import option_section
import numpy as np
from scipy.interpolate import UnivariateSpline

cosmo = names.cosmological_parameters


def setup(options):
    zmin = options[option_section, "zmin"]
    zmax = options[option_section, "zmax"]
    verbose = options[option_section, "verbose"]
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
    return verbose


def execute(block, config):

    verbose = config

    r, z, sigma_r = block.get_grid("sigma_r", "R", "z", "sigma2")
    omega_m = block[cosmo, "omega_m"]
    rho_c = 2.775e11

    def mass_given_r(r, omega_m):
        return 4.0 / 3.0 * np.pi * r**3 * omega_m * rho_c

    # z, R, blockname, crop_klim = config

    nz = len(z)
    nr = len(r)

    crit_r_array = np.zeros(nz)

    for iz, z_val in enumerate(range(nz)):
        sigma_r_at_z = sigma_r[:, iz]

        crit_density = 1.686**2
        sigma_r_reduced = sigma_r_at_z - crit_density
        f = UnivariateSpline(r, sigma_r_reduced, s=0)

        if verbose:
            print("This is redshift", z_val)
            print("Printing sigma_r_reduced")
            print(sigma_r_reduced)
            print("Printing roots of UnivariateSpline")
            print(f.roots())

        crit_r_array[iz] = f.roots()[0]

    crit_mass_array = mass_given_r(crit_r_array, omega_m)

    if verbose:
        print("Crit mass at z", list(zip(z, np.log10(crit_mass_array))))

    return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness
    return 0
