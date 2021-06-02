import os
from cosmosis.datablock import option_section, names
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import RectBivariateSpline
import scipy.special
import cluster_toolkit as ct
import matplotlib.pyplot as plt

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
        Rmis = options[section, 'rmis_array']

        return R_perp_min, R_perp_max, R_perp_bins, Radii_min, Radii_max, Radii_bins, M_min, M_max, M_bins, Rmis


def execute(block, config):
        R_perp_min, R_perp_max, R_perp_bins, Radii_min, Radii_max, Radii_bins, M_min, M_max, M_bins, Rmis= config
        
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

        mass_miscent_model = block['MassMiscentY1MortScalarIntegrand', "vals"]
        mass_cent_model = block['MassCentY1MortScalarIntegrand', "vals"]
        NC_miscent_model = block['NCRadMiscentY1MortScalarIntegrand', "vals"]
        NC_cent_model = block['NCCentY1MortScalarIntegrand', "vals"]
        miscent_fraction = block['cluster_abundance', 'fcen']
        masses = mass_miscent_model/NC_miscent_model
        truemass = mass_cent_model/NC_cent_model
        masses=masses.flatten()
        M_bins=len(masses+1)
        M=masses
        

        Xi_1 = np.zeros((M_bins, Radii_bins))
        Sigma_1 = np.zeros((M_bins, R_perp_bins))
        for i in range(M_bins):
                sigma1 = ct.deltasigma.Sigma_nfw_at_R(
                                R_perp, M[i], concentration, omega_m
                                )
                Xi_1[i] = ct.xi.xi_nfw_at_r(Radii, M[i], concentration, omega_m)
                Sigma_1[i] = sigma1
        zind=np.argmin((z-0.27)**2)
        Xi_2 = np.zeros(Radii_bins)
        ### define zind

        xi_mm = ct.xi.xi_mm_at_r(Radii, k_nl, P_k_nl[zind])
        Xi_2 = xi_mm
        Sigma_2 = ct.deltasigma.Sigma_at_R(
                                R_perp, Radii, xi_mm, 3.199267137797384375e+14, concentration, omega_m
                                )
        Bias = np.zeros(M_bins)
        for i in range(M_bins):
            Bias[i] = ct.bias.bias_at_M(M[i], k, P_k[zind], omega_m)

        logM = np.log(M)
        block["mass2sigma", "sigma_1"] = Sigma_1
        block["mass2sigma", "sigma_2"] = Sigma_2
        block["mass2sigma", "bias"] = Bias
        block["mass2sigma", "m_h"] = M
        block["mass2sigma", "z"] = z
        block["mass2sigma", "lnM"] = logM
        block["mass2sigma", "R_sigma_deltasigma"] = R_perp
        
        # define PR
        theta_array=np.linspace(0.0, np.pi * 2.0, 100)
        sigma_mis = np.zeros( [len(Rmis), len(R_perp)] )
        
        for ii, Rm in enumerate(Rmis):
           print(ii, Rm, np.log10(M[ii]), Bias[ii])
           print(Sigma_2)
           sigm_Rmis = np.zeros([len(theta_array), len(R_perp)])
           # set up fmis
           Sigma = Sigma_1[ii, :] + Bias[ii] * Sigma_2
           sigma_mis_single = ct.miscentering.Sigma_mis_single_at_R(R_perp, R_perp, Sigma, M[ii], concentration, omega_m, Rm) 
           sigma_mis[ii, :]=sigma_mis_single
       
        sigma_tot=np.zeros(len(R_perp))
        _tau=0.2
        P_Rerp= Rmis / _tau / _tau * np.exp(-Rmis / _tau)
        for jj in range(len(R_perp)):
            yy=sigma_mis[:, jj]* P_Rerp
            sigma_tot[jj] = np.trapz(yy, Rmis) 
        block["mass2sigma", "Averaged_sigma"] = sigma_tot
        plt.loglog(R_perp, sigma_tot, '--')
        plt.show()
        return 0




def cleanup(config):
        pass
