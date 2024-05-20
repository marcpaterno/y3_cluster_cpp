import os
from cosmosis.datablock import option_section, names
import numpy as np
from scipy.interpolate import interp1d, interp2d
from scipy.interpolate import RectBivariateSpline
import scipy.special
import cluster_toolkit as ct
from scipy import optimize


def setup(options):
    section = option_section
    return 0

def execute(block, config):
        test_mass=3.199267137797384375e+14
        test_z=0.2010101

        S2 = block["gammat", "Sigma_hh"] 
        S1 = block["gammat", "Sigma_nfw"]
        conc = block["gammat", "concentration"]
        Bias = block["gammat", "bias"]
    
        M = block["gammat", "m_h"]
        logM = block["gammat", "lnM"]
        zz = block["gammat", "z"]
        rhoc = block["gammat", "rhoc"]
        
        M = block["gammat", "m_h"]
        zz = block["gammat", "z"]
        logM = block["gammat", "lnM"]
        Rp = block["gammat", "r_sigma"]
        
        basepath = "%s/gammat" % os.environ["Y3_CLUSTER_CPP_DIR"]
        dat=np.genfromtxt('%s/test/Sigma_nfwonly.txt' % basepath, delimiter=',')
        test_r=dat[:, 0]; test_test=dat[:, 1]
        test_interp=interp2d(Rp, M, S1, kind='cubic')
        test_values=test_interp(test_r, test_mass)
        with open('%s/test/Sigma_nfwonly.out' % basepath, 'w') as outf:
            outf.write('r\t ytrue\t ytest\n')
            for ii in range(len(test_r)):
                outf.write('%f\t%.12e\t%.12e\n'%(test_r[ii], test_test[ii], test_values[ii]))

        dat=np.genfromtxt('%s/test/Sigma_mm_2halo.txt' % basepath, delimiter=',')
        test_r=dat[:, 0]; test_test=dat[:, 1]
        test_interp=interp2d(Rp, zz, S2, kind='quintic')
        test_values=test_interp(test_r, test_z)
        with open('%s/test/Sigma_mm_2halo.out' % basepath, 'w') as outf:
            outf.write('r\t ytrue\t ytest\n')
            for ii in range(len(test_r)):
                outf.write('%f\t%.12e\t%.12e\n'%(test_r[ii], test_test[ii], test_values[ii]))

        dat=np.genfromtxt('%s/test/bias.txt' % basepath, delimiter=',')
        test_mass=dat[:, 0]; test_test=dat[:, 2]
        test_z=dat[:, 1];
        test_interp=interp2d(zz, M, Bias, kind='cubic')
        test_value=test_interp(test_z, test_mass)
        with open('%s/test/bias.out' % basepath, 'w') as outf:
            outf.write('M\t z\t ytrue\t ytest\n')
            for ii in range(len(test_mass)):
                test_value=test_interp(test_z[ii], test_mass[ii])
                outf.write('%f\t%f\t%.12e\t%.12e\n'%(test_mass[ii], test_z[ii], test_test[ii], test_value))

        return 0

def check(x,label):
    print('Check Variable: %s'%label,x)
    print('Shape of vector %s is: '%label, x.shape)

def cleanup(config):
    pass
