import os
import numpy as np
from cosmosis.datablock import option_section, names
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import interp1d
#import matplotlib.pyplot as plt

def setup(options):
     cons = 'test'
     return cons

def execute(block, config):
    """
    Tinker bias

    Convention:
        mass is in units of M_sun
        m200 is in units of 10^14 M_sun
    """
    cons = config
    cosmo = names.cosmological_parameters
    dis = names.distances
    gro = "growth_parameters"
    tbf = "tinker_bias_function"

    #cosmosis_dir = "output_cross/"
    #z_table = np.genfromtxt(os.path.join(cosmosis_dir,"growth_parameters/z.txt"), skip_header=1)
    #growth_table = np.genfromtxt(os.path.join(cosmosis_dir,"growth_parameters/d_z.txt"), skip_header=1)
    z_table = block[gro,'z']
    growth_table = block[gro,'d_z']
    growth_function = interp1d(z_table, growth_table)

    #cosmosis_dir = "output_cross/"
    #z_table = np.genfromtxt(os.path.join(cosmosis_dir,"sigma_r/z.txt"), skip_header=1)
    #r_table = np.genfromtxt(os.path.join(cosmosis_dir,"sigma_r/r.txt"), skip_header=1)
    #sigmar_table=np.genfromtxt(os.path.join(cosmosis_dir,"sigma_r/sigma2.txt"), skip_header=1)
    z_table = block["sigma_r", "z"]
    r_table = block["sigma_r", "r"]
    sigmar_table=block["sigma_r","sigma2"]
    sigmar_table=np.sqrt(sigmar_table)

    # tinker mass function is defined as  M_critical; 
    # the arithmetric below results in M_critical
    #       to get M_background mult by omega_m
    omega_m = block[cosmo,"omega_m"]
    h0 = block[cosmo,"h0"]
    omega_l = 1. - omega_m
    rho_m=2.775e11 # h^2 M_solar Mpc^-3.

    m_table=(4.0*3.1416/3.0)*rho_m*r_table**3  # M_solar/h
    #2D interpolator into mass function; mass in linear in solar masses
    ln_m_table = np.log(m_table)
    sigma_mass = RectBivariateSpline(ln_m_table,z_table,sigmar_table)

    y = np.log10(200.)
    A = 1.0 + 0.24 * y * np.exp(-(4/y)**4)
    alpha = 0.44 * y - 0.88
    B = 0.183
    beta = 1.5
    C = 0.019 + 0.107*y + 0.19*np.exp(-(4/y)**4)
    ceta = 2.4
    delta = 1.686

    # let's evaulate it for the places we need
    bias_table = np.array([])
    ln_m_interp_h = np.log(10.0**(np.linspace(7, 17, 100))/h0)
    for zed in z_table:
        # tinker mass function is defined as  M_critical; 
        # the arithmetric below results in M_critical
        #       to get M_background mult by omega_m
        E0 =  (omega_l + omega_m * (1. + zed)**3)
        omega_mz = omega_m*(1.0 + zed)**3/E0 #omega_m(z)
        m_table=(4.0*3.1416/3.0)*rho_m*omega_mz*r_table**3  # M_solar/h
        ln_m_table = np.log(m_table)

        sigma_mass = RectBivariateSpline(ln_m_table,z_table,sigmar_table)  
        nu = delta/sigma_mass(ln_m_interp_h,zed)

        bias = 1 - A*(nu**alpha)/(nu**alpha + delta**alpha) + \
            B* nu**beta + C* nu**ceta
        bias = bias.flatten()
        if bias_table.size == 0 :
            bias_table = bias
        else :
            bias_table = np.vstack([bias_table, bias])

    block[tbf,'z'] = z_table
    block[tbf,'ln_mass_h'] = ln_m_interp_h
    block[tbf,'bias'] = bias_table

    return 0
def cleanup(config):
    #nothing to do here!  We just include this 
    # for completeness
    return 0

