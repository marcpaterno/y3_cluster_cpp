import os
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section, names
import numpy as np
from scipy.interpolate import interp1d, interp2d
from scipy.interpolate import RectBivariateSpline
import scipy.special
import cluster_toolkit as ct
from scipy import optimize
section_name="y1_analysis_mor"
import matplotlib.pyplot as plt


def setup(options):
    section = option_section
    cov=np.genfromtxt('y1_data/Cov_ij_bestfit_DESY1_105.txt')
    NCdata=np.genfromtxt('y1_data/NC_DESY1A1_REDMAPPER_3zbin_4lbin_Mcut1.e13_LM_v6.4.17_MISCCorr.dat')
    Mdata=np.genfromtxt('y1_data/LogWLMass_4lbin_3zbin_redmapper_desy1_v3.dat')


    #row1:  z_1 lambda_1
    #row2:  z_1 lambda_2
    #row3:  z_1 lambda_3
    Om_corr = np.genfromtxt('y1_data/logMdOm_DESY1.dat')
    lnAS_corr = np.genfromtxt('y1_data/logMdlnAs_DESY1.dat')
    return cov, NCdata, Mdata, Om_corr, lnAS_corr



def execute(block, config):

        covmat, ncdata, Mdata, Om_corr, lnAS_corr = config

        nn = (block[section_name, "n_vals"])
        logM = np.log10(block[section_name, "totm_vals"]/nn)
        Omega_m = (block["cosmological_parameters", "omega_m"])
        loge10As = (block["cosmological_parameters", "log1e10As"])

        l_low = [20.0, 30.0, 45.0, 60.0]
        l_high = [30.0, 45.0, 60.0, 190.0]
        z_low = [0.2, 0.35, 0.5]
        z_high = [0.35, 0.5, 0.65]
       
        theory_nc=np.array([])
        theory_M=np.array([])
        corr=np.array([])
        for jj in range(len(z_low)):
            for ii in range(len(l_low)):
                theory_nc=np.append(theory_nc, nn[ii, jj])#/(l_high[ii]-l_low[ii]))
                theory_M=np.append(theory_M, logM[ii, jj])
                roww=jj*4+ii
                corr_jjii=Om_corr[roww]*(Omega_m-0.3)+lnAS_corr[roww]*(loge10As-2.98239371)
                corr=np.append(corr, corr_jjii)
        richnesses=ncdata[:, 2]
        redshifts=ncdata[:, 0]
        nc_data=ncdata[:, 4]
        M_data=Mdata[:, 4]
        #M_data=M_data+corr

        data=np.append(nc_data, M_data)
        theory=np.append(theory_nc, theory_M)
        delta = data - theory
        weight = np.linalg.inv(covmat)
        loglike1 = -0.5 * np.dot(delta, np.dot(weight, delta))
        print("theory vector M and NC:", theory_M, theory_nc)
        print("Data vector M and Corr:", M_data, corr)
  
        '''
        delta2 = mass_data - theory_M
        weight2 = np.linalg.inv(covmat2)
        loglike2 = -0.5 * np.dot(delta2, np.dot(weight2, delta2))
        '''

        block[section_names.likelihoods, 'Y1AnalysisLike_like'] = loglike1
        #plt.plot(richnesses, nc_data, 'ko')
        #plt.plot(richnesses, theory_nc, 'ro')
        #plt.show()

        return 0


def cleanup(config):
	pass
