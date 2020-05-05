import os
from cosmosis.datablock import option_section, names
import numpy as np
from scipy.interpolate import interp1d
import astropy.io.fits as pyfits

cosmo = names.cosmological_parameters

def setup(options):
	section = option_section
        z_file_name=os.path.expandvars('${Y3_CLUSTER_CPP_DIR}/data/data_DESY1/test_cluster_Y1.fits'
        z_data=pyfits.open(z_file_name)[6].data
        z_array=z_data['z_mid']
        prob_array=z_data['bin1']

        z_min = options[section,"z_min"]
        z_max = options[section,"z_max"]
        z_bins = int(options[section,"z_bins"])
        zl_array=np.linspace(z_min, z_max, z_bins)

	return z_array, prob_array, zl_array

def execute(block, config):
        z_vals, prob_vals, zl_array=config

        G=4.51710305e-48 #Mpc^3 / M_sol / s^2
        c=9.71561189e-15 #Mpc/s
	h0 = block[cosmo, "h0"]
	delta_z = block["photoz", "delta_z"]
        pzs=interp1d(z_vals, prob_vals)
        print np.sum(prob_vals)

        z_da=block['distances', 'z']
        da=block['distances', 'd_a'] # Mpc
        da_func=interp1d(z_da, da)
        da_zl=da_func(zl_array)
        da_zs=da_func(z_vals)

        int_sci=np.zeros(len(zl_array))
        for jj, zt in enumerate(zl_array):
            da_zl_zs = da_zs - (1.0+zt)/(1.0+z_vals) *da_zl[jj]
            int_array= np.maximum(np.zeros(len(z_vals)), pzs(z_vals+delta_z)*4.0*np.pi*G/c/c * da_zl[jj] * da_zl_zs / da_zs)
            
            int_sci[jj]=np.trapz(int_array, z_vals)*h0
        block["average_sigma_crit_inv", "zlense"]=zl_array
        block["average_sigma_crit_inv", "sci_average"]=int_sci
        return 0


def cleanup(config):
	pass
