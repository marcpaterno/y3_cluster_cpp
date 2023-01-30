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
# - The data vector is (NC, GT, Wp)
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
# Wp stands for cluster correlation function (3d or 2d)
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

names = ['NC', 'Shear', 'Wp']

import h5py
def readHDF(fname, path):
    master = h5py.File(fname,'r')
    return master[path][:][:]

def readDataVector(fname, path='data'):
    mydict = dict().fromkeynames(names)
    for name in names:
        mydict[name] = readHDF(fname, name+'/'+path)
    return mydict

def setup(options):
    section = option_section
    datavector_fname = 0

    dataDict = readDataVector(datavector_fname, 'data')
    invCovDict = readDataVector(datavector_fname, 'invcov')
    return dataDict, invCovDict


def execute(block, config):
    # read from the data
    dataDict, invCovDict datavector_fname = config
    
    # pull predictions from the datablock; 
    # arrays with shape (Nrbins x Nlbdins x Nzbins, )
    # [(r0, l0, z0), (r1, l0, z0), (r2, l0, z0), ..., (r0, zn, lm), (r1, zn, lm), ...]
    # For NC; Nrbins = 0 
    NC  = block["numberCounts", "vals"].flatten()
    GT  = block["shear", "shear_cen"].flatten()
    WP = block["wp", "wp_cen"].flatten()
    # Sum Operator: take the average by dividing Nc[1]

    theoryDict = {'NC': NC, 'Shear': GT, 'Wp': WP}

    # logLikelihood
    # block-diagonal covariances with R
    logLike = 0
    for name in names:
        data = dataDict[name]
        theory = theoryDict[name]

        # the current covariance is block-diagonal
        # i.e. no covariance between redshift and lambda bins
        invcov = invCovDict[name]

        delta = data - theory
        logLike += -0.5 * np.dot(delta, np.dot(invcov, delta))

    # put into the datablock
    block[section_names.likelihoods, 'log_like'] = loglike
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
