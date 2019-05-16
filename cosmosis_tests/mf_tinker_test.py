import os
from cosmosis.datablock import option_section, names
import numpy as np
from scipy.interpolate import interp1d, interp2d
from scipy.interpolate import RectBivariateSpline
import scipy.special
import cluster_toolkit as ct
from scipy import optimize
section_name="mass_function"

def setup(options):
    section = option_section
    return 0


def execute(block, config):

        Omega_M = block["cosmological_parameters", "omega_M"]-block["cosmological_parameters", "omega_nu"]
        masses = np.log10(Omega_M*block[section_name, "m_h"])
        zz = block[section_name, "z"]
	nm = (block[section_name, "dndlnmh"])

        dat=np.genfromtxt('/cosmosis/cosmosis-standard-library/y3_cluster_cpp/test/test_hmf_z0_z03.dat')
        test_z= dat[:, 0] 
        test_m= dat[:, 1] 
        test_test=dat[:, 2]*(10.0**test_m)

        test_interp=interp2d(masses, zz, nm, kind='quintic')
        with open('/cosmosis/cosmosis-standard-library/y3_cluster_cpp/data/masses_test.out', 'w') as outf:
            outf.write('z\t mass\t ytrue\t ytest\n')
            for ii in range(len(test_m)):
                test_value = test_interp(test_m[ii], test_z[ii])*(4.50732047e-02* (test_m[ii] - 13.8124426028) + 1.01958078e+00)
                outf.write('%f\t%f\t%e\t%e\n'%(test_z[ii], test_m[ii], test_test[ii], test_value))


        return 0


def cleanup(config):
	pass
