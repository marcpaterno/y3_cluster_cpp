import os
from cosmosis.datablock import option_section, names
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import RectBivariateSpline
import scipy.special
import cluster_toolkit as ct
import massconcen

cosmo = names.cosmological_parameters


def setup(options):
    section = option_section

    #saveOutput = int(options[section, "saveOutput"])

    #Mpc/h comoving distance, distance on the sky
    R_perp_min = options[section,"R_perp_min"]
    R_perp_max = options[section,"R_perp_max"]
    R_perp_bins = int(options[section,"R_perp_bins"])

    # radii of xi_hm, in Mpc/h
    Radii_min = options[section,"Radii_min"]
    Radii_max = options[section,"Radii_max"]
    Radii_bins = int(options[section,"Radii_bins"])

    #mass (float): Halo mass Msun/h.
    M_min = options[section,"M_min"]
    M_max = options[section,"M_max"]
    M_bins = int(options[section,"M_bins"])

    return R_perp_min, R_perp_max, R_perp_bins, Radii_min, Radii_max, Radii_bins, M_min, M_max, M_bins


def execute(block, config):
    R_perp_min, R_perp_max, R_perp_bins, Radii_min, Radii_max, Radii_bins, M_min, M_max, M_bins= config

    omega_m = block[cosmo, "omega_m"]
    omega_b = block[cosmo, "omega_b"]
    sigma8 = block[cosmo, "sigma_8"]
    h0 = block[cosmo, "h0"]

    mstar = block["mstar", "mstar"]
    mstar_z = block["mstar", "z"]

    k = block["matter_power_lin", "k_h"]
    P_k = block["matter_power_lin", "p_k"]
    z_k = block["matter_power_lin", "z"]

    k_nl = block["matter_power_nl", "k_h"]
    P_k_nl = block["matter_power_nl", "p_k"]

    z = block["matter_power_nl", "z"]
    nz = len(z)

    R_perp = np.logspace(np.log10(R_perp_min), np.log10(R_perp_max), R_perp_bins)
    Radii = np.logspace(np.log10(Radii_min), np.log10(Radii_max), Radii_bins)
    M = np.logspace(np.log10(M_min), np.log10(M_max), M_bins)

    Xi_1 = np.zeros((M_bins, Radii_bins))
    Sigma_1 = np.zeros((M_bins, R_perp_bins))
    Fx = np.zeros((M_bins, R_perp_bins))
    for i in range(M_bins):
        # the mass concentration relation is a function of redshift, but the way Y. Zhang designed 
        # this, the 1-halo term is independent of redshift, which is computationally efficient.
        # We'll have to assume a mean redshift and compute the c from M at that z.
        zfix = 0.33
        concentration = massconcen.c_from_m200 (M[i], zfix, omega_m, omega_b, sigma8, h0, mstar, mstar_z)
        sigma1 = ct.deltasigma.Sigma_nfw_at_R(
            R_perp, M[i], concentration, omega_m)
        Xi_1[i] = ct.xi.xi_nfw_at_r(Radii, M[i], concentration, omega_m)
        Sigma_1[i] = sigma1

        # we'll need the F(x) of the projected NFW for the boost factor. sigma1 contains it
        # but also has a normalizatin factor that we don't want. (sigma1 is in Mpc/h)
        # What we'll do is to reproduce the norm in ct.C_deltasigma.c's Sigma_nfw_at_R_arr
        # at lines 52-55.
        # Note that the gx computation has a singularity at x=1, and will blow up there. 
        # WE SHOULD FIX THAT.
        delta = 200 # M200m
        rhocrit = 2.77533742639e+11  #h^2 M_sun Mpc^-3
        # 1e4*3.*Mpcperkm*Mpcperkm/(8.*PI*G); units are SM h^2/Mpc^3
        rhom = omega_m*rhocrit # SM h^2/Mpc^3
        deltac = delta*(1./3.)*concentration**3/ (
            np.log(1+concentration)-concentration/(1+concentration) )
        Rdelta = (M[i]/(1.3333*np.pi*rhom*delta))**(1./3.)
        Rscale = Rdelta/concentration
        norm = 2*Rscale*deltac*rhom*1e-12
        Fx[i] = sigma1/norm

    # 2-halo term = FT(PS)*bias, but here we don't calculate the bias
	# concentration is a dummy parameter in this call, passed all the way down to c-code
	# Sigma_at_R_Arr (in C_deltasigma.c) , put into a params structure at line 161, 
	# and not used in the integration.
    dummy_concentration = 5.
    Xi_2 = np.zeros((nz, Radii_bins))
    Sigma_2 = np.zeros((nz, R_perp_bins))
    for i in range(nz):
        xi_mm = ct.xi.xi_mm_at_r(Radii, k_nl, P_k_nl[i])
        Xi_2[i] = xi_mm
        sigma2 = ct.deltasigma.Sigma_at_R(
            R_perp, Radii, xi_mm, 3.2e+14, dummy_concentration, omega_m)
        Sigma_2[i] = sigma2

    # calculate the bias
    # FIXME make bias real
    Bias = np.zeros((M_bins, nz))
    for i in range(M_bins):
        for j in range(nz):
             Bias[i][j] = ct.bias.bias_at_M(M[i], k, P_k[j], omega_m)

    logM = np.log(M)
    block["deltasigma", "Xi_1"] = Xi_1
    block["deltasigma", "Xi_2"] = Xi_2
    block["deltasigma", "R_Xi"] = Radii
    block["deltasigma", "sigma_1"] = Sigma_1
    block["deltasigma", "sigma_2"] = Sigma_2
    block["deltasigma", "bias"] = Bias
    block["deltasigma", "m_h"] = M
    block["deltasigma", "z"] = z
    block["deltasigma", "lnM"] = logM
    block["deltasigma", "R_sigma_deltasigma"] = R_perp
    block["deltasigma", "Fx"] = Fx
    return 0


def cleanup(config):
    pass
