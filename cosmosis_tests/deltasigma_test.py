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

        Xi_1 = block["deltasigma", "Xi_1"]
        Xi_2 = block["deltasigma", "Xi_2"]
        Radii = block["deltasigma", "R_Xi"]
        S1 = block["deltasigma", "sigma_1"]
        S2 = block["deltasigma", "sigma_2"]
        Bias = block["deltasigma", "bias"]
        M = block["deltasigma", "m_h"]
        zz = block["deltasigma", "z"]
        logM = block["deltasigma", "lnM"]
        RDS = block["deltasigma", "R_sigma_deltasigma"]
        test_mass=3.199267137797384375e+14
        test_z=0.2010101

        basepath = "%s/deltasigma" % os.environ["Y3_CLUSTER_CPP_DIR"]
        dat=np.genfromtxt('%s/test/xi_nfwonly.txt' % basepath, delimiter=',')
        test_r=dat[:, 0]; test_test=dat[:, 1]
        test_interp=interp2d(Radii, M, Xi_1, kind='cubic')
        test_values=test_interp(test_r, test_mass)
        with open('%s/test/xi_nfwonly.out' % basepath, 'w') as outf:
            outf.write('r\t ytrue\t ytest\n')
            for ii in range(len(test_r)):
                outf.write('%f\t%f\t%f\n'%(test_r[ii], test_test[ii], test_values[ii]))

        dat=np.genfromtxt('%s/test/Sigma_nfwonly.txt' % basepath, delimiter=',')
        test_r=dat[:, 0]; test_test=dat[:, 1]
        test_interp=interp2d(RDS, M, S1, kind='cubic')
        test_values=test_interp(test_r, test_mass)
        with open('%s/test/Sigma_nfwonly.out' % basepath, 'w') as outf:
            outf.write('r\t ytrue\t ytest\n')
            for ii in range(len(test_r)):
                outf.write('%f\t%f\t%f\n'%(test_r[ii], test_test[ii], test_values[ii]))


        dat=np.genfromtxt('%s/test/xi_mm_2halo.txt' % basepath, delimiter=',')
        test_r=dat[:, 0]; test_test=dat[:, 1]
        test_interp=interp2d(Radii, zz, Xi_2, kind='cubic')
        test_values=test_interp(test_r, test_z)
        with open('%s/test/xi_mm_2halo.out' % basepath, 'w') as outf:
            outf.write('r\t ytrue\t ytest\n')
            for ii in range(len(test_r)):
                outf.write('%f\t%f\t%f\n'%(test_r[ii], test_test[ii], test_values[ii]))

        dat=np.genfromtxt('%s/test/Sigma_mm_2halo.txt' % basepath, delimiter=',')
        test_r=dat[:, 0]; test_test=dat[:, 1]
        test_interp=interp2d(RDS, zz, S2, kind='quintic')
        test_values=test_interp(test_r, test_z)
        with open('%s/test/Sigma_mm_2halo.out' % basepath, 'w') as outf:
            outf.write('r\t ytrue\t ytest\n')
            for ii in range(len(test_r)):
                outf.write('%f\t%f\t%f\n'%(test_r[ii], test_test[ii], test_values[ii]))

        dat=np.genfromtxt('%s/test/bias.txt' % basepath, delimiter=',')
        test_mass=dat[:, 0]; test_test=dat[:, 2]
        test_z=dat[:, 1];
        test_interp=interp2d(zz, M, Bias, kind='cubic')
        test_value=test_interp(test_z, test_mass)
        with open('%s/test/bias.out' % basepath, 'w') as outf:
            outf.write('M\t z\t ytrue\t ytest\n')
            for ii in range(len(test_mass)):
                test_value=test_interp(test_z[ii], test_mass[ii])
                outf.write('%f\t%f\t%f\t%f\n'%(test_mass[ii], test_z[ii], test_test[ii], test_value))


        return 0


def cleanup(config):
    pass
