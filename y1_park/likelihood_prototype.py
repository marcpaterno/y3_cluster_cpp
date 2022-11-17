import os
from cosmosis.datablock import option_section, names
import numpy as np
from scipy.special import erf
import cluster_toolkit as ct

###########################################################
############ Likelihood: Prototype #############
# First likelihood version of the DES Y1 Park model project
# In this version we implement:
# - The data vector is (NC, GT, CRF)
# - The covariance matrix is computed by Jackknife
# - We use mock data build from the Buzzard Sims v1.9.8
###########################################################
############ Codes To Validate ###############
# There are new scripts we should validate
# 1) We build kappa instead of Sigma (see code.cc)
# 2) We compute shear from the radial kappa excess (see buildShear.py)
# 4) We use Sigma_crit computed from the effective source redshift (see buildShear.py)
# 5) tbd
# Note: We don't model boost-factor. 
###########################################################
###########################################################
# Useful Nomeclature
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

    # centering fraction
    fcen = (block["cluster_abundance", "fcen"])
    
    ## load model predictions
    # centered component
    NC_cen = block["NCCentY1MortScalarIntegrand", "vals"]    
    GT_cen  = block["MassCentY1MortScalarIntegrand", "vals"]
    CRF_cen= block["CorrCentY1MortScalarIntegrand", "vals"]
    
    # mis-centered component
    NC_mis = block["NCMiscentY1MortScalarIntegrand", "vals"]
    GT_mis  = block["MassMiscentY1MortScalarIntegrand", "vals"]
    CRF_mis= block["CorrMiscentY1MortScalarIntegrand", "vals"]
    
    # to check results; debug is a global variable
    if debug:
        print('number count predictions - centered')
        print(NC_cen)
        
        print('number count predictions - mis-centered')
        print(NC_mis)
    
    # cosmological parameters
    Omega_m = (block["cosmological_parameters", "omega_m"])
    loge10As = (block["cosmological_parameters", "log1e10As"])
    fcen = (block["cluster_abundance", "fcen"])
    
    # final predictions
    NC  = fcen*NC_cen + (1.0-fcen)*NC_mis
    GT  = fcen*GT_cen + (1.0-fcen)*GT_mis
    CRF = fcen*CRF_cen+ (1.0-fcen)*CRF_mis
    
    # why?????
    logM = np.log10(GT/NC)
    
    # put into the datablock
    block["final_model", "NC"] = NC
    block["final_model", "GT"] = GT
    block["final_model", "CRF"] = CRF
    return 0

def cleanup(config):
    pass
