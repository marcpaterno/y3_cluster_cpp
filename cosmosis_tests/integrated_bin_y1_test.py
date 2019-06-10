import os
from cosmosis.datablock import option_section, names
import numpy as np
from scipy.interpolate import interp1d, interp2d
from scipy.interpolate import RectBivariateSpline
import scipy.special
import cluster_toolkit as ct
from scipy import optimize
section_name="y1_analysis"

def setup(options):
    section = option_section
    return 0


def execute(block, config):

        nn = (block[section_name, "n_vals"][:, 0])
        logM = np.log10(block[section_name, "totm_vals"][:, 0]/nn)

        l_low = [20.0, 27.9, 37.6, 50.3]
        l_high = [27.9, 37.6, 50.3, 69.3]
        teo_nc = [3698.17826207, 1780.16346057, 923.24586279, 468.39002918]
        teo_logm = [14.17440261, 14.35812886, 14.52711919, 14.69816962]
        with open('/cosmosis/cosmosis-standard-library/y3_cluster_cpp/data/integrated_bin_y1_test.out', 'w') as outf:
            outf.write('llo\t lhi \t yntrue\t yntest \t ymtrue\t ymtest\n')
            for ii in range(len(l_low)):
                outf.write('%f\t%f\t%f\t%f\t%f\t%f\n'%(l_low[ii], l_high[ii], teo_nc[ii], nn[ii], teo_logm[ii], logM[ii]))


        return 0


def cleanup(config):
	pass
