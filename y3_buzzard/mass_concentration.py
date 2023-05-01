from __future__ import print_function
from cosmosis.datablock import names
from cosmosis.datablock import option_section
import numpy as np
from scipy.interpolate import UnivariateSpline
import numpy as np
import hydro_mc
from scipy.interpolate import interp1d


# m200 is mass M200m in linear units of solar mass/h^-1
def c_from_m200 (m200, z, omega_m, omega_b, sigma8, h0, mstar,mstar_z) :
    # using parameters from individual, all row of table 1
    A =  3.44
    b =  430.49
    m = -0.10
    co =  3.19
    # row 3, stacked NFW
    #A =  4.61
    #b =  638.65
    #m = -0.07
    #co =  3.59
    # row 4, stacked einasto
    #A =  63.2
    #b =  431.48
    #m = -0.01
    #co =  3.36
    # now M_* should be calculated using equations 13-17, but we will approximate it
    # reproduces Colossus' Child18 model reasonably well
    #Mstarc =  10**(14.76 * .808**z)

    # from Child's eq 13 via mstar.py module
    interp=interp1d(mstar_z, mstar)
    Mstar = interp(z)
    #z = 0.01 # c = 5.02 at 14.0
    #z = 0.33 # c = 4.60 at 14.0
    #z = 0.66 # c = 4.20 at 14.0
    #z = 0.99 # c = 3.90 at 14.0

    # eq 18 demands M200c, not M200m. So we will use the conversion equations in
    # Ragagnin, Saro, Singh, & Dolag 2021 and code at 
    # https://github.com/aragagnin/hydro_mc/blob/master/examples/sample_mm.py
    # so we will need to import hydro_mc
    # https://github.com/aragagnin/hydro_mc/blob/master/hydro_mc.py
    # Then:
    M_200c = hydro_mc.mass_from_mm_relation(
        '200m','200c', M=m200, a=1./(1+z),
        omega_m= omega_m, omega_b= omega_b, sigma8= sigma8, h0= h0)
    # note M_200c should be strictly less then M_200m

    # Now we do eq 18
    mmb = M_200c/Mstar/b
    c200 = A * ( mmb**m * (1+mmb)**(-1*m) -1) + co
    #print(np.log10(m200),np.log10(Mstar),m200/Mstar,c200)

    return c200


# use the mass concentration relation from
    # Ragagnin, Saro, Singh, & Dolag 2021 and code at 
    # so we will need to import hydro_mc
    # https://github.com/aragagnin/hydro_mc/blob/master/hydro_mc.py
    #
# m200 is mass M200m in linear units of solar mass/h^-1
def c_from_m200_ragagnin (m200, z, omega_m, omega_b, sigma8, h0) :

    z= 0.5; m200=0.3e15
    cs = c_from_m200_child (m200, z, omega_m, omega_b, sigma8, h0) 

    c200 = hydro_mc.concentration_from_mc_relation(
            '200m', M=m200, a=1./(1+z), omega_m= omega_m, omega_b= omega_b, sigma8= sigma8, h0= h0)

    print("===================",c200, cs)
    import sys
    sys.tracebacklimit = 0
    raise Exception("here")
    return c200



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

    R, z, sigma_r = block.get_grid("sigma_r", "R", "z", "sigma_r")
    print(z)
    omega_m = block[cosmo, "omega_m"]
    rho_c = 2.775e11

    print(sigma_r)

    def mass_given_r(r, omega_m):
        return 4.0 / 3.0 * np.pi * r**3 * omega_m * rho_c

    # z, R, blockname, crop_klim = config

    nz = len(z)
    nr = len(R)

    print(nr, nz)
    print(sigma_r.shape)

    crit_r_array = np.zeros(nz)

    for iz, z_val in enumerate(z):

        sigma_r_at_z = sigma_r[:, iz]

        crit_density = 1.686
        sigma_r_reduced = sigma_r_at_z - crit_density
        f = UnivariateSpline(R, sigma_r_reduced, s=0)

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
