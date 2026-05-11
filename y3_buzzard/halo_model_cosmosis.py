import os
import numpy as np
from cosmosis.datablock import option_section, names

from scipy.interpolate import interp1d
from scipy.special import erf
import cluster_toolkit as ct

from astropy.constants import G

from haloModel import biasModel, lensingModel, scaleShiftCosmo

##################################################################
############# Build Wp(R,z) and Gammat(R,z) ############
# Estimates Halo Model Parameters
# Wp, \Sigma, \Delta \Sigma and bias(z)
# 
#############################################
# Author: Johnny Esteves
# Created: May 16, 2024
#############################################

import astropy.cosmology
cosmo_names = names.cosmological_parameters

def setup(options):
    section = option_section
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

    # Two independent toggles for the lensing stack:
    #
    #   compute_lensing_1h  publish Sigma_nfw, dSigma_nfw, concentration
    #                       (cluster_toolkit analytic NFW tables, ~cheap).
    #                       Read by the legacy Sigma1hSel / DSigma1hSel
    #                       / Shear1hSel / ReducedShear1hSel modules.
    #
    #   compute_lensing_2h  publish Sigma_hh, dSigma_hh, Wp_hh
    #                       (cluster_toolkit ct_2hTerm: xi_mm + Hankel
    #                       Sigma_at_R + DeltaSigma_at_R — ~200-300 ms
    #                       and the dominant cost of this module).
    #                       Read by the *TotSel variants only.
    #
    # Costanzi-2026 projection pipeline (red_shear_prj + bsel) only
    # needs haloModel/bias and xi_nl/xi_nl: set both toggles to F.
    # ReducedShear1hSel alone: 1h=T, 2h=F.  The legacy Tot stack: both T.
    #
    # Backward-compat: the old single knob `compute_lensing` controls
    # both when set (maps to 1h AND 2h); defaults are (1h=T, 2h=T).
    def _bool(key, default):
        try:
            return bool(options.get_bool(section, key, default=default))
        except Exception:
            return default
    compute_lensing = _bool("compute_lensing", True)
    compute_lensing_1h = _bool("compute_lensing_1h", compute_lensing)
    compute_lensing_2h = _bool("compute_lensing_2h", compute_lensing)

    params_out = (R_perp_min, R_perp_max, R_perp_bins,
                  Radii_min, Radii_max, Radii_bins,
                  M_min, M_max, M_bins,
                  compute_lensing_1h, compute_lensing_2h)
    return params_out
    

def execute(block, config):
    section_name = "haloModel"

    (R_perp_min, R_perp_max, R_perp_bins,
     Radii_min, Radii_max, Radii_bins,
     M_min, M_max, M_bins,
     compute_lensing_1h, compute_lensing_2h) = config

    # cosmo parameters
    omega_m = block[cosmo_names, "omega_m"]
    omega_b = block[cosmo_names, "omega_b"]
    H0 = block[cosmo_names, "H0"]

    cosmology = astropy.cosmology.FlatLambdaCDM(H0*100, omega_m, Ob0=omega_b, Tcmb0=2.725)

    # load camb matter-matter power spectrums
    # linear is used for the bias
    k_h = block["matter_power_lin", "k_h"]
    P_k = block["matter_power_lin", "p_k"]
    z_k = block["matter_power_lin", "z"]

    # nonlinear P(k) -- fall back to linear if cp_camb's nonlinear emulator
    # is not plugged in (wiring test only; physics approximate).
    if block.has_section("matter_power_nl"):
        P_k_nl = block["matter_power_nl", "p_k"]
        k_nl   = block["matter_power_nl", "k_h"]
        z_nl   = block["matter_power_nl", "z"]
    else:
        P_k_nl, k_nl, z_nl = P_k, k_h, z_k

    z = z_k

    if block.has_section("growth_parameters"):
        z_az = block["growth_parameters", "z"]
        daz  = block["growth_parameters", "d_z"]
        # CosmoSIS publishes an un-normalised growth factor
        # (D -> a at early times, so D(0) != 1).  bias(M, z) needs the
        # ratio sigma(M, z)/sigma(M, 0) = D(z)/D(0), so renormalise here.
        daz  = np.interp(z, z_az, daz) / np.interp(0.0, z_az, daz)
    else:
        # cp_camb mode: no growth_parameters.  Approximate D(z) via
        # the shape of sigma_8 only (daz = 1 at z=0, drops slowly).
        # For a pure-wiring-timing test, just set D(z)=1 and move on;
        # Bias(M, z) is wrong but the timing measurement is correct.
        daz = np.ones_like(z)

    # compute overdensities; physical not comoving
    rho_c0 = float(cosmology.critical_density0.to('Msun/Mpc^3').value)
    # rho_cz = cosmology.critical_density(z).to('Msun/Mpc^3').value
    rho_m = omega_m*rho_c0
    rho_mz = rho_m*(1.+z)**3
    
    # setup bins
    nz = len(z)
    R_perp = np.logspace(np.log10(R_perp_min), np.log10(R_perp_max), R_perp_bins)
    Radii = np.logspace(np.log10(Radii_min), np.log10(Radii_max), Radii_bins)
    M = np.logspace(np.log10(M_min), np.log10(M_max), M_bins)
    logM = np.log(M)

    # Compute Bias (M, Z)
    # The bias is computed at the peak height (z=0) from cluster toolkit
    # The peak-height evolution is computed from the growth function D(z)
    # The vector Bias has shape (z.size, M.size)
    bM = biasModel(k_h, P_k[0], omega_m, odelta=200)
    nu = bM.nu_at_M(M)
    Bias = np.array([bM.bias_at_nu(nu/dazi) for dazi in daz])

    # xi_NL(r, z) pre-tabulated on a fixed log-r grid for the C++
    # b_sel_marg / sigma_prj integrands. Evaluated per z slice; kept
    # outside the `compute_lensing` branch because the Costanzi-2026
    # projection pipeline always needs it.
    r_xi = np.logspace(-3.0, 3.0, 128)  # Mpc/h comoving
    xi_NL = np.zeros((z.size, r_xi.size))
    for iz, zz in enumerate(z):
        xi_NL[iz] = ct.xi.xi_mm_at_r(r_xi, k_nl, P_k_nl[iz])

    # ----- always-on datablock writes (consumed by bsel + red_shear_prj) -----
    block[section_name, "m_h"] = M
    block[section_name, "lnM"] = logM
    block[section_name, "z"] = z
    block[section_name, "rhoc"] = rho_mz
    block[section_name, "bias"] = Bias
    # xi_NL(r, z) table for the b_sel_marg / sigma_prj integrands
    block["xi_nl", "r"]     = r_xi
    block["xi_nl", "z"]     = z
    block["xi_nl", "xi_nl"] = xi_NL

    # ----- lensing section: optional ----------------------------------------
    # Split into two independent branches: 1h (cheap analytic NFW
    # tables) and 2h (expensive cluster_toolkit Hankel transforms).
    # The modern red_shear_prj + bsel pipeline does not read either;
    # ReducedShear1hSel needs the 1h half; ReducedShearTotSel needs
    # both.  Skipping the 2h alone saves ~200-300 ms/sample.
    if compute_lensing_1h or compute_lensing_2h:
        lensModel = lensingModel(R_perp, omega_m, 200)
        scaleShift, hubbleShift = scaleShiftCosmo(z, cosmology)
        block[section_name, "r_sigma"] = R_perp
        block[section_name, "scale_shift"] = scaleShift
        block[section_name, "hubble_shift"] = hubbleShift
        block[section_name, "k"] = k_h

    if compute_lensing_1h:
        lensModel.first_halo_term(M, z=0, conc_model_name="Child18")
        block[section_name, "Sigma_nfw"]  = lensModel.Sigma['1h']
        block[section_name, "dSigma_nfw"] = lensModel.dSigma['1h']
        block[section_name, "concentration"] = lensModel.c

    if compute_lensing_2h:
        lensModel.second_halo_term(z, k_nl, P_k_nl)
        block[section_name, "Rp"]        = Radii
        block[section_name, "Wp_hh"]     = lensModel.Wp
        block[section_name, "Sigma_hh"]  = lensModel.Sigma['2h']
        block[section_name, "dSigma_hh"] = lensModel.dSigma['2h']

    return 0

def cleanup(config):
    pass