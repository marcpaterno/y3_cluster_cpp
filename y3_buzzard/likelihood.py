import os
from cosmosis.datablock import option_section, names
import numpy as np

###########################################################
################### UNDER CONSTRUCTION #################
# ToDos:
#   - create read data function
#   - check data structure of the theory outpu

###########################################################
################### Likelihood: Prototype #################
# First likelihood version of the DES Y1 Park model project
#
# Likelihood = (dataVector - model) cov_inv (dataVector - model)^T
#
# In this version we implement:
# - The data vector is (NC, GT, CRF)
# - The covariance matrix is computed by Jackknife
# - We use mock data build from the Buzzard Sims v1.9.8
#
###########################################################
################### Codes To Validate #####################
# There are new scripts we should validate
# 1) We build kappa instead of Sigma (see code.cc)
# 2) We compute shear from the radial kappa excess (see buildShear.py)
# 4) We use Sigma_crit computed from the effective source redshift (see buildShear.py)
# 5) tbd
# Note: We don't model boost-factor. 
###########################################################
###########################################################
##### Useful Nomeclature
###########################################################
# Datavector codenames
# ---------------------------------------------------------
# NC stands for cluster number counts 
# GT stands for gamma_t (shear)
# CRF stands for cluster correlation function (3d or 2d)
# ---------------------------------------------------------
# Miscentering model codenames
# _cen for centered models
# _mis for mis-centered models
# fcen is the centered fraction
# ---------------------------------------------------------
###########################################################
# Author: Johnny Esteves
# Created: 17th Nov, 2022
###########################################################

cosmo = names.cosmological_parameters
debug = True

def setup(options):
    section = option_section

    return 0

def execute(block, config):
    # script based on y1_rerun/y1_analysis_mor_wmiscent_likelihood.py

    # read from the data
    covmat, ncdata, Mdata, Om_corr, lnAS_corr = config
    
    # pull predictions from the datablock; 
    # arrays with shape (Nlbdins, Nzbins)
    NC  = block["summedNumbersCentBu", "vals"]
    GT  = block["avg_shear", "shear_cen"]
    CRF = block["avgPkDampCentBu", "vals"]

    # lambda and redshift bins
    # why are they not setup in a global variable?
    l_low = [20.0, 30.0, 45.0, 60.0]
    l_high = [30.0, 45.0, 60.0, 190.0]
    z_low = [0.2, 0.35, 0.5]
    z_high = [0.35, 0.5, 0.65]

    # concatanates the datavector over lambda and redshift bins
    # NC has shape of (N_lbins, N_zbins)
    # put models as 1d array
    theory_nc = concatenate(NC, lbins=len(l_low), zbins=len(z_low))
    theory_gt = concatenate(GT, lbins=len(l_low), zbins=len(z_low))
    theory_crf = concatenate(CRF, lbins=len(l_low), zbins=len(z_low))
    
    ## check this lines on the original code
    # corr_jjii=Om_corr[roww]*(Omega_m-0.3)+lnAS_corr[roww]*(loge10As-2.98239371)
    # corr=np.append(corr, corr_jjii)
    # M_data=M_data+corr

    # unfolds the data
    # data vector 
    data_nc = ncdata[:, 4]
    # dummy variables; to check data structrue
    data_gt = gtdata[:, 4]
    data_crf = crfdata[:, 4]

    # joint datavector
    data = np.append(data_nc, data_gt)
    data = np.append(data, data_crf)

    # model predictions
    theory = np.append(theory_nc, theory_gt)
    theory = np.append(theory, theory_crf)

    # covariance matrix
    cov_inv =np.linalg.inv(covmat)

    # logLikelihood
    delta = data - theory
    loglike = -0.5 * np.dot(delta, np.dot(cov_inv, delta))

    # put into the datablock
    block[section_names.likelihoods, 'BuzY3'] = loglike
    return 0

def concatenate_datavector(x, lbins=4, zbins=3):
    # concatanates the datavector over lambda and redshift bins
    # x has shape of (N_lbins, N_zbins)
    # returns 1d array 
    y = np.array([])
    for jj in range(len(N_zbins)):
        for ii in range(len(N_lbins)):
            y = np.append(y, x[ii, jj])
    return y
    
def cleanup(config):
    pass
