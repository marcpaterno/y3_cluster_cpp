import os
import numpy as np
from scipy.interpolate import interp1d
from cosmosis.datablock import option_section, names
from scipy import integrate

# redshift mean with format [z0, z1, z2, z0, z1, z2, ...]
from setup_bins import zmeans_ij

from haloModel import biasModel
from nfwModel import deltaSigmaNFW_Analytical
# from cluster_toolkit import averaging


#############################################
#############################################
# Author: Johnny Esteves
# Created: 17th Nov, 2022
#############################################

cosmo_names = names.cosmological_parameters

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
    method = str(options[section, "method"])
    lbdbins = np.array(options[section, "lbdbins"])
    rbins = np.array(options[section, "radii"])
    zlow = np.array(options[section, "zo_low"])
    zhigh = np.array(options[section, "zo_high"])
    zmean = 0.5*(zlow+zhigh)

    lbdMeans = lbdbins[:-1]+np.diff(lbdbins)/2.
    meffParams = []#TBD 
    ceffParams = []#TBD
    return method, rbins, zmean, lbdMeans, meffParams, ceffParams


def execute(block, config):
    """Compute an effective lensing profile

    Grab the 2h term from haloModel
    Compute the 1h term for Meff
    Compute the bias for Meff and zeff
    """
    method, rbins, zeff, lbdMeans, meffParams, ceffParams = config
    zvec = block["haloModel", "z"]
    rhomz = block["haloModel", "rhoc"]
    omega_m = block[cosmo_names, "omega_m"]
    rhom0 = rhomz[0]

    z_az = block["growth_parameters", "z"]
    _daz = block["growth_parameters", "d_z"]
    # daz = np.interp(z, z_az, _daz)

    k_h = block["matter_power_lin", "k_h"]
    P_k = block["matter_power_lin", "p_k"]

    Rp = block["haloModel","r_sigma"]
    dSigma2h = block["haloModel", "dSigma_hh"] 

    sigma_crit_inv_r = block["sigmaCritInv", "sigma_crit_inv"]
    # take only the values that does not depend on boost-factor
    sigma_crit_inv = sigma_crit_inv_r[-1]

    # grab mass vector values
    Meff = np.logspace(11,12.5, 12)
    params = [1.,1.,1.]
    ceff = powerLaw(Meff, 0., params)
    bM = biasModel(k_h, P_k[0], omega_m, odelta=200)
    nuEff = bM.nu_at_M(Meff)

    rhomz_eff = np.interp(zeff, zvec, rhomz)
    da_eff = np.interp(zeff, z_az, _daz)
    sigma_crit_inv_eff = np.interp(zeff, zvec, sigma_crit_inv)

    lensEff = []
    meff = []
    ceff =[]
    iz = np.interp(zeff, zvec, np.arange(zvec.size, dtype=np.int32)).astype(np.int32)
    dSigma2h_eff = dSigma2h[iz]
    for lbd in lbdMeans:
        M = powerLaw(lbd, zeff, meffParams)
        nu = bM.nu_at_M(M)
        bias = bM.bias_at_nu(nu/da_eff)
        c = powerLaw(M, zeff, ceffParams)
        dsigma1h = deltaSigmaNFW_Analytical(rbins, M, c, rho_c=rhomz_eff)
        dS = hmModel(dsigma1h, dSigma2h_eff, bias)
        shear = dS*sigma_crit_inv_eff
        lensEff.append(shear)
        meff.append(M)
        ceff.append(c)

    # put into the datablock
    block["opSelBias", "r"] = rbins 
    block["opSelBias", "ceff"] = ceff
    block["opSelBias", "meff"] = meff
    block["opSelBias", "shear"] = np.array(lensEff)

    return 0

def get_params(block,name="meff"):
    a = block[f"{name}-alpha"]
    b = block[f"{name}-beta"]
    c = block[f"{name}-gamma"]
    return [a,b,c]

def hmModel(halo1, halo2, bias):
    return np.max(halo1, bias*halo2, axis=0)

def powerLaw(x, z, params):
    a, b, c, x0, z0 = params
    return np.log(a) + b*np.log(x/x0) + c*np.log((1+z)/(1+z0))

def cleanup(config):
    pass
