import os
import numpy as np

#####################################################################
###################### Setup Bin Scheme  ############################
# Defines three redshift bins and 4 lambda bins
#
# Uses this set of global variables along the code
###################################################################
# Author: Johnny Esteves
# Created: 18th Nov, 2022
###################################################################

# redshift bins
# zbins = np.array([0.2, 0.35, 0.5, 0.65])
# zmeans = np.array([0.275, 0.4, 0.575])
# zmeans_ij = np.array([0.275, 0.4, 0.575, 0.275, 0.4, 0.575, 0.275, 0.4, 0.575, 0.275, 0.4, 0.575])

# Buzzard Mock
zbins = np.array([0.2, 0.32, 0.373, 0.51, 0.65])
zmin_list = np.array([0.2, 0.373, 0.51])
zmax_list = np.array([0.32, 0.51, 0.64])
zmeans = np.array([0.25, 0.44, 0.575])
zmeans_ij = np.array([0.25, 0.44, 0.575, 0.25, 0.44, 0.575, 0.25, 0.44, 0.575, 0.25, 0.44, 0.575])

# lambda bins
lbdbins =  np.array([5, 15, 25, 40, 160])
lbdmeans = 0.5*(lbdbins[1:]+lbdbins[:-1])
lbdmeans_ij = np.repeat(lbdmeans, 3)


