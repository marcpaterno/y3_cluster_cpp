import os
from cosmosis.datablock import option_section, names
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import RectBivariateSpline
import scipy.special
import cluster_toolkit as ct
import matplotlib.pyplot as plt
from scipy.integrate import simps

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
        #M_min = options[section,"M_min"]
        #M_max = options[section,"M_max"]
        #M_bins = int(options[section,"M_bins"])
        Rmis_array = options[section, 'rmis_array']
        zo_low = options[section,"zo_low"]
        zo_high = options[section,"zo_high"]

        return R_perp_min, R_perp_max, R_perp_bins, Radii_min, Radii_max, Radii_bins, Rmis_array, zo_low, zo_high


def execute(block, config):
        R_perp_min, R_perp_max, R_perp_bins, Radii_min, Radii_max, Radii_bins, Rmis_array, zo_low, zo_high= config
        
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

        mass_miscent_model = block['MassRadMiscentY1MortScalarIntegrand', "vals"]
        mass_cent_model = block['MassCentY1MortScalarIntegrand', "vals"]
        NC_miscent_model = block['NCRadMiscentY1MortScalarIntegrand', "vals"]
        NC_cent_model = block['NCCentY1MortScalarIntegrand', "vals"]
        cent_fraction = block['cluster_abundance', 'fcen']
        masses = mass_miscent_model/NC_miscent_model
        truemasses = mass_cent_model/NC_cent_model
        #masses=masses.flatten()
        print (masses.shape, NC_miscent_model.shape, len(masses))

        ### define zind
        zl, zl_indices = np.unique(zo_low, return_index=True)
        zh, zh_indices = np.unique(zo_high, return_index=True)
        Rmis, Rmis_indices = np.unique(Rmis_array, return_index=True)
        print(zl_indices, zh_indices)
        if not (np.array_equal(zl_indices, zh_indices) & (len(zl)*len(Rmis) == len(Rmis_array) ) ):
            raise Exception('check if zo_low, zo_high, Rmis pairs are matched properly')
        
        
        sigma_mis = np.zeros( [len(masses)*len(zl_indices), len(R_perp)] )
        sigma_cent = np.zeros( [len(masses)*len(zl_indices), len(R_perp)] )
        sigma_tot = np.zeros( [len(masses)*len(zl_indices), len(R_perp)] )
        DeltaSigma = np.zeros( [len(masses)*len(zl_indices), len(R_perp)] )
        
        masses_tot = np.zeros([len(masses), len(zl_indices)])
        masses_miscent = np.zeros([len(masses), len(zl_indices)])
        masses_miscent_tot = np.zeros([len(masses), len(zl_indices)])
        masses_cent = np.zeros([len(masses), len(zl_indices)])
        ncs_tot = np.zeros([len(masses), len(zl_indices)])
        ncs_miscent = np.zeros([len(masses), len(zl_indices)])
        ncs_cent = np.zeros([len(masses), len(zl_indices)])
        for ww in range(len(zl_indices)):
            zmid=(zl[ww]+zh[ww])/2.0
            zind=np.argmin((z-zmid)**2)
            xi_mm = ct.xi.xi_mm_at_r(Radii, k_nl, P_k_nl[zind])
            ww_indmin=ww*len(Rmis)
            ww_indmax=(ww+1)*len(Rmis)
         
            for jj in range(len(masses)):
                sigma_mis_arr = np.zeros( [len(Rmis), len(R_perp)] )
                for ii, Rm in enumerate(Rmis):
                    #print(ii, Rm, np.log10(M[ii]), Bias[ii])
                    #Sigma = Sigma_1[ii, :] + Bias[ii] * Sigma_2
                    #sigma_mis_single = ct.miscentering.Sigma_mis_single_at_R(R_perp, R_perp, Sigma, M[ii], concentration, omega_m, Rm) 
                    Xi_1 = ct.xi.xi_nfw_at_r(Radii, masses[jj, ww_indmin+ii], concentration, omega_m)
                    Bias = ct.bias.bias_at_M(masses[jj, ww_indmin+ii], k, P_k[zind], omega_m)
                    xi_tot = ct.xi.xi_hm(Xi_1, Bias * xi_mm)
                    #xi_tot =  Bias * xi_mm
                    Sigma = ct.deltasigma.Sigma_at_R(R_perp, Radii, xi_tot, 3.199267137797384375e+14, concentration, omega_m)
                    sigma_mis_single = ct.miscentering.Sigma_mis_single_at_R(R_perp, R_perp, Sigma, masses[jj, ww_indmin+ii], concentration, omega_m, Rm) 
                    sigma_mis_arr[ii, :]=sigma_mis_single
       
                #sigma_mis=np.zeros(len(R_perp))
                _tau=0.2
                P_Rerp= Rmis / _tau / _tau * np.exp(-Rmis / _tau)
                for kk in range(len(R_perp)):
                    yy=sigma_mis_arr[:, kk]* P_Rerp
                    sigma_mis[jj*len(zl_indices)+ww, kk] = np.trapz(yy, Rmis)/np.trapz(P_Rerp, Rmis)

                masses_miscent_tot[jj, ww] = np.trapz(masses[jj, ww_indmin:ww_indmax]*NC_miscent_model[jj, ww_indmin:ww_indmax] *P_Rerp, Rmis)/np.trapz(P_Rerp, Rmis) 
                masses_cent[jj, ww] = truemasses[jj, ww_indmin] 
                ncs_miscent[jj, ww] = np.trapz(NC_miscent_model[jj, ww_indmin:ww_indmax] * P_Rerp, Rmis)/np.trapz(P_Rerp, Rmis) 
                ncs_cent[jj, ww] = NC_cent_model[jj, ww_indmin] 
                ncs_tot[jj, ww] = cent_fraction * NC_cent_model[jj, ww_indmin] + (1.0-cent_fraction)*ncs_miscent[jj, ww] 
                masses_miscent[jj, ww] = masses_miscent_tot[jj, ww]/ncs_miscent[jj, ww] 
                masses_tot[jj, ww] = ( cent_fraction*(truemasses[jj, ww_indmin]*NC_cent_model[jj, ww_indmin]) + (1.0-cent_fraction)*masses_miscent_tot[jj, ww])/ncs_tot[jj, ww] 
            

                #centered clusters
                Xi_1 = ct.xi.xi_nfw_at_r(Radii, truemasses[jj, ww_indmin], concentration, omega_m)
                Bias = ct.bias.bias_at_M(truemasses[jj, ww_indmin], k, P_k[zind], omega_m)
                xi_tot = ct.xi.xi_hm(Xi_1, Bias * xi_mm)
                sigma_cent[jj*len(zl_indices)+ww] = ct.deltasigma.Sigma_at_R(R_perp, Radii, xi_tot, 3.199267137797384375e+14, concentration, omega_m)
                #total averaged of miscentered, and centered
                sigma_tot[jj*len(zl_indices)+ww] = cent_fraction * sigma_cent[jj*len(zl_indices)+ww] + (1.0-cent_fraction)*sigma_mis[jj*len(zl_indices)+ww]
                DeltaSigma[jj*len(zl_indices)+ww] = ct.deltasigma.DeltaSigma_at_R(R_perp, R_perp, sigma_tot[jj*len(zl_indices)+ww], masses_tot[jj, ww], concentration, omega_m)

                ''' 
                plt.loglog(R_perp, sigma_tot[jj*len(zl_indices)+ww], 'k')
                plt.loglog(R_perp, sigma_mis[jj*len(zl_indices)+ww], 'b:')
                plt.loglog(R_perp, sigma_cent[jj*len(zl_indices)+ww], 'r:')
          
                test_sigma=np.genfromtxt('sigma_y1_output/sigmacenty1mortscalarintegrand/vals.txt')            
                test_nc=np.genfromtxt('sigma_y1_output/nccenty1mortscalarintegrand/vals.txt')            
                rad = [0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]
                print(np.shape(test_nc[jj]), np.shape(test_sigma[jj]), np.shape(rad))
                plt.loglog(rad, test_sigma[jj]/test_nc[jj], 'g--')
                test_sigma=np.genfromtxt('sigma_y1_output/sigmamiscenty1mortscalarintegrand/vals.txt') 
                test_nc=np.genfromtxt('sigma_y1_output/ncmiscenty1mortscalarintegrand/vals.txt')            
                plt.loglog(rad, test_sigma[jj]/test_nc[jj], 'g:')
                '''
            
        #block["mass2sigma", "sigma_1"] = Sigma_1
        #block["mass2sigma", "sigma_2"] = Sigma_2
        #block["mass2sigma", "bias"] = Bias
        block["mass2sigma", "NC_cent"] = ncs_cent
        block["mass2sigma", "NC_miscent"] = ncs_miscent
        block["mass2sigma", "NC_tot"] = ncs_tot
        block["mass2sigma", "mcent_h"] = masses_cent # * NC_cent_model
        block["mass2sigma", "m_miscent_h"] = masses_miscent # * ncs_miscent
        block["mass2sigma", "m_average_h"] = masses_tot
        block["mass2sigma", "z"] = z
        block["mass2sigma", "R_sigma_deltasigma"] = R_perp
        block["mass2sigma", "Averaged_sigma"] = sigma_tot
        block["mass2sigma", "delta_sigma"] = DeltaSigma
        block["mass2sigma", "Sigma_cent"] = sigma_cent
        block["mass2sigma", "Sigma_miscent"] = sigma_mis
        plt.show()
        return 0




def cleanup(config):
        pass
