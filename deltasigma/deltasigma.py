import os
from cosmosis.datablock import option_section, names
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import RectBivariateSpline
import scipy.special
import cluster_toolkit as ct

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
        concentration = block["cluster_mass_profile", "concentration"]
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
        for i in range(M_bins):
                sigma1 = ct.deltasigma.Sigma_nfw_at_R(
                                R_perp, M[i], concentration, omega_m
                                )
                Xi_1[i] = ct.xi.xi_nfw_at_r(Radii, M[i], concentration, omega_m)
                Sigma_1[i] = sigma1

        Xi_2 = np.zeros((nz, Radii_bins))
        Sigma_2 = np.zeros((nz, R_perp_bins))
        for i in range(nz):
                xi_mm = ct.xi.xi_mm_at_r(Radii, k_nl, P_k_nl[i])
                Xi_2[i] = xi_mm
                sigma2 = ct.deltasigma.Sigma_at_R(
                                R_perp, Radii, xi_mm, 3.199267137797384375e+14, concentration, omega_m
                                )
                Sigma_2[i] = sigma2

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
        return 0


def cleanup(config):
        pass
