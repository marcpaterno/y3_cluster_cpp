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
	DS1 = block["deltasigma", "deltasigma_1"]
	DS2 = block["deltasigma", "deltasigma_2"]
	S1 = block["deltasigma", "sigma_1"]
	S2 = block["deltasigma", "sigma_2"]
	Bias = block["deltasigma", "bias"]
	M = block["deltasigma", "m_h"]
	zz = block["deltasigma", "z"]
	logM = block["deltasigma", "lnM"]
	RDS = block["deltasigma", "R_sigma_deltasigma"]
        test_mass=10.0**14.5
        test_z=0.201010

        dat=np.genfromtxt('test/xi_14.5_nfw.txt', delimiter=',')
        test_r=dat[:, 0]; test_test=dat[:, 1]
        test_interp=interp2d(Radii, M, Xi_1, kind='cubic')
        test_values=test_interp(test_r, test_mass)
        with open('test/xi_14.5_nfw.out', 'w') as outf:
            outf.write('r\t ytrue\t ytest\n')
            for ii in range(len(test_r)):
                outf.write('%f\t%f\t%f\n'%(test_r[ii], test_test[ii], test_values[ii]))

        dat=np.genfromtxt('test/Sigma_14.5_nfwonly.txt', delimiter=',')
        test_r=dat[:, 0]; test_test=dat[:, 1]
        test_interp=interp2d(RDS, M, S1, kind='cubic')
        test_values=test_interp(test_r, test_mass)
        with open('test/Sigma_14.5_nfw.out', 'w') as outf:
            outf.write('r\t ytrue\t ytest\n')
            for ii in range(len(test_r)):
                outf.write('%f\t%f\t%f\n'%(test_r[ii], test_test[ii], test_values[ii]))

        dat=np.genfromtxt('test/DS_14.5_nfw.txt', delimiter=',')
        test_r=dat[:, 0]; test_test=dat[:, 1]
        test_interp=interp2d(RDS, M, DS1, kind='cubic')
        test_values=test_interp(test_r, test_mass)
        with open('test/DS_14.5_nfw.out', 'w') as outf:
            outf.write('r\t ytrue\t ytest\n')
            for ii in range(len(test_r)):
                outf.write('%f\t%f\t%f\n'%(test_r[ii], test_test[ii], test_values[ii]))

        dat=np.genfromtxt('test/xi_2halo.txt', delimiter=',')
        test_r=dat[:, 0]; test_test=dat[:, 1]
        test_interp=interp2d(Radii, zz, Xi_2, kind='cubic')
        test_values=test_interp(test_r, test_z)
        with open('test/xi_2halo.out', 'w') as outf:
            outf.write('r\t ytrue\t ytest\n')
            for ii in range(len(test_r)):
                outf.write('%f\t%f\t%f\n'%(test_r[ii], test_test[ii], test_values[ii]))

        dat=np.genfromtxt('test/Sigma_14.5_2halo.txt', delimiter=',')
        test_r=dat[:, 0]; test_test=dat[:, 1]
        test_interp=interp2d(RDS, zz, S2, kind='cubic')
        test_values=test_interp(test_r, test_z)
        with open('test/Sigma_14.5_2halo.out', 'w') as outf:
            outf.write('r\t ytrue\t ytest\n')
            for ii in range(len(test_r)):
                outf.write('%f\t%f\t%f\n'%(test_r[ii], test_test[ii], test_values[ii]))

        dat=np.genfromtxt('test/DS_14.5_2halo.txt', delimiter=',')
        test_r=dat[:, 0]; test_test=dat[:, 1]
        test_interp=interp2d(RDS, zz, DS2, kind='cubic')
        print RDS.shape, zz.shape, DS2.shape
        test_values=test_interp(test_r, test_z)
        with open('test/DS_14.5_2halo.out', 'w') as outf:
            outf.write('r\t ytrue\t ytest\n')
            for ii in range(len(test_r)):
                outf.write('%f\t%f\t%f\n'%(test_r[ii], test_test[ii], test_values[ii]))


        dat=np.genfromtxt('test/bias.txt', delimiter=',')
        test_mass=dat[:, 0]; test_test=dat[:, 2]
        test_z=dat[:, 1];
        print zz.shape, M.shape, Bias.shape
        test_interp=interp2d(zz, M, Bias, kind='cubic')
        test_value=test_interp(test_z, test_mass)
        with open('test/bias.out', 'w') as outf:
            outf.write('M\t z\t ytrue\t ytest\n')
            for ii in range(len(test_mass)):
                test_value=test_interp(test_z[ii], test_mass[ii])
                outf.write('%f\t%f\t%f\t%f\n'%(test_mass[ii], test_z[ii], test_test[ii], test_value))


        return 0


def cleanup(config):
	pass
