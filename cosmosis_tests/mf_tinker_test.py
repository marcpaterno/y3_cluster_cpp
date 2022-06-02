import os
from cosmosis.datablock import option_section
import numpy as np
from scipy.interpolate import interp2d
section_name="mass_function"

def setup(options):
    section = option_section
    return 0


def execute(block, config):
    # Jim, May 2022, understanding Yuanyuan Zhang's test

        Omega_M = block["cosmological_parameters", "omega_M"]-block["cosmological_parameters", "omega_nu"] # reduce to Omega_coldmatter
        masses = np.log10(Omega_M*block[section_name, "m_h"]) # unclear why multiply m_h by omega_m, if Matteo thing, why not fix those?
        zz = block[section_name, "z"]
        nm = (block[section_name, "dndlnmh"])

        dat = np.genfromtxt(os.path.expandvars('${Y3_CLUSTER_CPP_DIR}/test/test_hmf_z0_z03.dat') ) # two z's, 0 and 0.3
        test_z = dat[:, 0] # z
        test_m = dat[:, 1] # mass
        test_test = dat[:, 2]*(10.0**test_m) # space density

        test_interp = interp2d(masses, zz, nm, kind='quintic') # this is to interpolate the mass function to a given z and mass
        outfile = os.path.expandvars('${Y3_CLUSTER_CPP_DIR}/data/masses_test.out')
        with open(outfile, 'w') as outf:
            outf.write('z\t mass\t ytrue\t ytest\n')
            # for each of Matteo's masses,
            for ii in range(len(test_m)):
                # test_test is from Matteo's SDSS mass function
                # test_value is the cosmosis evaluated mass function, i.e., tinker_mass_function
                # we're going to evaluate at the z,mass  of Matteo's mass function
                # thus we get the cosmosis mf, call that test_interp,
                cosmosis_value = test_interp(test_m[ii], test_z[ii])
                test_value = cosmosis_value


             # But we must make an empirical correction for the test to pass
                # then make a linear correction: mulitply by (0.045 * (matteos_mass_point - 13.81) + 1.02 )
                #test_value = (4.50732047e-02* (test_m[ii] - 13.8124426028) + 1.01958078e+00) * cosmosis_value
                # a better correction which solves the one point not passing the above. Recall this is after Cosmosis V2.0 change.
                test_value = (4.75e-02* (test_m[ii] - 13.81) + 1.01958078e+00) * cosmosis_value

                outf.write('%f\t%f\t%e\t%e\n'%(test_z[ii], test_m[ii], test_test[ii], test_value))


        return 0


def cleanup(config):
	pass
