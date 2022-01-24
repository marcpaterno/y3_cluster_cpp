import os
from cosmosis.datablock import option_section, names
import numpy as np
from scipy.interpolate import interp1d, interp2d
from scipy.interpolate import RectBivariateSpline
import scipy.special
import cluster_toolkit as ct
from scipy import optimize
section_name="y1_analysis"
import matplotlib.pyplot as plt

def setup(options):
    section = option_section
    return 0


def execute(block, config):

        nn = (block[section_name, "n_vals"])
        logM = np.log10(block[section_name, "totm_vals"]/nn)

        l_low = [20.0, 30.0, 45.0, 60.0]
        l_high = [30.0, 45.0, 60.0, 190.0]
        z_low = [0.2, 0.35, 0.5]
        z_high = [0.35, 0.5, 0.65]
        teo_nc = np.asarray([[1304.9894721052153, 528.95563650496740, 155.52915961522930, 96.457909451135279],
                  [1905.3087824221111, 699.03200511844079, 187.99368993235004, 107.53455525604019],
                  [1919.6887249119727,  643.53225321401715, 156.62336704176681, 80.214466253146725]])

        teo_logm = np.asarray([[14.066856527866555, 14.296261606649280, 14.487104492397755, 14.717555086743570],
            [14.055662163683669, 14.282850922869070, 14.470959036325867, 14.682116585781188],
            [14.027073034542163, 14.254233304065675, 14.443080091917629, 14.641826179274396]])
        outfile = os.path.expandvars("${Y3_CLUSTER_CPP_DIR}/data/integrated_bin_y1_test.out")
        with open(outfile, 'w') as outf:
            outf.write('llo\t lhi\t zlo\t zhi\t yntrue\t yntest \t ymtrue\t ymtest\n')
            for ii in range(len(l_low)):
               for jj in range(len(z_low)):
                outf.write('%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n'%(l_low[ii], l_high[ii], z_low[jj], z_high[jj], teo_nc[jj, ii], nn[ii, jj]/(l_high[ii]-l_low[ii]), teo_logm[jj, ii], logM[ii, jj]))
        dat=np.genfromtxt(outfile)
        plt.plot(dat[:, 0], dat[:, 4], 'ko')
        plt.plot(dat[:, 0], dat[:, 5], 'ro')
        plt.show()
        plt.plot(dat[:, 0], dat[:, 6], 'ks')
        plt.plot(dat[:, 0], dat[:, 7], 'rs')
        plt.show()
        return 0


def cleanup(config):
	pass
