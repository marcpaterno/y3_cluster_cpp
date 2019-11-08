import os
from cosmosis.datablock import option_section
import numpy as np
from scipy.interpolate import interp2d
section_name="mass_function"

def setup(options):
    section = option_section
    return 0


def execute(block, config):

        Omega_M = block["cosmological_parameters", "omega_M"]-block["cosmological_parameters", "omega_nu"]
        masses = np.log10(Omega_M*block[section_name, "m_h"])
        zz = block[section_name, "z"]
        nm = (block[section_name, "dndlnmh"])

        dat = np.genfromtxt(os.path.expandvars('${COSMOSIS_SRC_DIR}/cosmosis-standard-library/y3_cluster_cpp/test/test_hmf_z0_z03.dat') )
        test_z = dat[:, 0] 
        test_m = dat[:, 1] 
        test_test = dat[:, 2]*(10.0**test_m)

        test_interp = interp2d(masses, zz, nm, kind='quintic')
        outfile = os.path.expandvars('${COSMOSIS_SRC_DIR}/cosmosis-standard-library/y3_cluster_cpp/data/masses_test.out')
        with open(outfile, 'w') as outf:
            outf.write('z\t mass\t ytrue\t ytest\n')
            for ii in range(len(test_m)):
                test_value = test_interp(test_m[ii], test_z[ii])*(4.50732047e-02* (test_m[ii] - 13.8124426028) + 1.01958078e+00)
                outf.write('%f\t%f\t%e\t%e\n'%(test_z[ii], test_m[ii], test_test[ii], test_value))


        return 0


def cleanup(config):
	pass
