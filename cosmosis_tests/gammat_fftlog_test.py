import os
from cosmosis.datablock import option_section
import numpy as np
from scipy.interpolate import interp2d
section_name="mass_function"

def setup(options):
    section = option_section
    ztest  = block[section, "ztest"]
    return ztest


def execute(block, config):
    # load data
    ztest = config
    z     = block["correlationFunction", "z"]

    # get redshift slice
    # ixz = np.interp(zt, z, np.arange(z.size,dtype=int)).astype(int)

    # deltaSigma Profiles output
    mh  = block["correlationFunction", "m_h"]
    lnM = block["correlationFunction", "lnM"]
    rperp = block["correlationFunction", "r_sigma"]

    sigma1 = block["correlationFunction", "Sigma_nfw"] # radii, mass
    sigma2 = block["correlationFunction", "Sigma_hh"]  # radii, redshift
    
    # 2d interpolation
    test_interp = interp2d(lnM, rperp, sigma1, kind='quintic')

    # get benchmark data
    dat = np.genfromtxt(os.path.expandvars('${Y3_CLUSTER_CPP_DIR}/test/test_hmf_z0_z03.dat') ) # two z's, 0 and 0.3
    test_rperp = dat[:, 0] # rperp
    test_mass = dat[:, 1] # mass
    test_z    = dat[:, 2] # redshift
    test_test = dat[:, 2] # nfw-profile

    # write output file
    outfile = os.path.expandvars('${Y3_CLUSTER_CPP_DIR}/data/sigma1_test.out')
    with open(outfile, 'w') as outf:
        outf.write('rperp\t mass\t ytrue\t ytest\n')
        # for each of Matteo's masses,
        for ii in range(len(test_m)):
            # test_test is from Matteo's SDSS mass function
            # test_value is the cosmosis evaluated mass function, i.e., tinker_mass_function
            # we're going to evaluate at the z,mass  of Matteo's mass function
            # thus we get the cosmosis mf, call that test_interp,
            cosmosis_value = test_interp(np.log(test_m[ii])), test_rperp)
            test_value = cosmosis_value
            outf.write('%f\t%f\t%e\t%e\n'%(test_rperp[ii], test_m[ii], test_test[ii], test_value))
        return 0


def cleanup(config):
	pass
