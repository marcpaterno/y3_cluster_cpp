import os
import numpy as np
from cosmosis.datablock import option_section, names
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import interp1d

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
    tbf = "tinker_bias_function"

    z_table = block["sigma_r", "z"]
    r_table = block["sigma_r", "r"]
    sigmar_table=block["sigma_r","sigma2"]
    sigmar_table=np.sqrt(sigmar_table)
    omega_m = block[cosmo,"omega_m"]
    rho_c=2.7754e11 # h^2 M_solar Mpc^-3.

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
    nu_table = np.array([])
    ln_m_interp_h = np.log( 10.0**(np.linspace(7, 17, 100)) )
    for zed in z_table:
        # tinker bias function is defined as  M_mean; https://arxiv.org/pdf/1001.3162.pdf 
        omega_mz = omega_m*(1.0 + zed)**3
        m_table=(4.0*3.1416/3.0)*rho_c*omega_mz*r_table**3  # M_solar/h
        ln_m_table = np.log(m_table)

        sigma_mass = RectBivariateSpline(ln_m_table,z_table,sigmar_table)  
        nu = delta/sigma_mass(ln_m_interp_h,zed)

        bias = 1 - A*(nu**alpha)/(nu**alpha + delta**alpha) + \
            B* nu**beta + C* nu**ceta
        bias = bias.flatten()
        nu = nu.flatten()
        if bias_table.size == 0 :
            bias_table = bias
            nu_table = nu
        else :
            bias_table = np.vstack([bias_table, bias])
            nu_table = np.vstack([nu_table, nu])

    block[tbf,'z'] = z_table
    block[tbf,'ln_mass_h'] = ln_m_interp_h
    block[tbf,'bias'] = bias_table
    block[tbf,'nu'] = nu_table

    return 0
def cleanup(config):
    #nothing to do here!  We just include this 
    # for completeness
    return 0

