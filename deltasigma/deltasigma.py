import os
from cosmosis.datablock import option_section, names
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import RectBivariateSpline
import scipy.special
import cluster_toolkit as ct

cosmo = names.cosmological_parameters

def setup(options):
	section = option_section

	#saveOutput = int(options[section, "saveOutput"])

	#Mpc/h comoving distance, distance on the sky
	R_perp_min = options[section,"R_perp_min"]
	R_perp_max = options[section,"R_perp_max"]
	R_perp_bins = int(options[section,"R_perp_bins"])

	# radii of xi_hm, in Mpc/h
	Radii_min = options[section,"Radii_min"]
	Radii_max = options[section,"Radii_max"]
	Radii_bins = int(options[section,"Radii_bins"])

	#mass (float): Halo mass Msun/h.
	M_min = options[section,"M_min"]
	M_max = options[section,"M_max"]
	M_bins = int(options[section,"M_bins"])

	concentration = options[section,"concentration"]

	return  R_perp_min, R_perp_max, R_perp_bins, Radii_min, Radii_max, Radii_bins, M_min, M_max, M_bins, concentration

def execute(block,config):
	R_perp_min, R_perp_max, R_perp_bins, Radii_min, Radii_max, Radii_bins, M_min, M_max, M_bins, concentration = config

	omega_m = block[cosmo,"omega_m"]
	n_s = block[cosmo,"n_s"]
	omega_b = block[cosmo,"omega_b"]
	h0 = block[cosmo,"h0"]

	k = block["matter_power_lin", "k_h"]
	P_k = block["matter_power_lin", "p_k"]

	k_nl = block["matter_power_nl", "k_h"]
	P_k_nl = block["matter_power_nl", "p_k"]

	z= block["matter_power_lin", "z"]
	nz=len(z)

	R_perp = np.logspace(np.log10(R_perp_min),np.log10(R_perp_max),R_perp_bins)
	Radii = np.logspace(np.log10(Radii_min),np.log10(Radii_max),Radii_bins)
	M = np.logspace(np.log10(M_min),np.log10(M_max),M_bins)

	Size1 = M_bins * R_perp_bins
	Deltasigma_1 = np.ndarray(Size1)
	Deltasigma_1.shape = (M_bins,R_perp_bins)

	for i in range(M_bins):
		Sigma1 = ct.deltasigma.Sigma_nfw_at_R(R_perp, M[i], concentration, omega_m)
		Deltasigma1 = ct.deltasigma.DeltaSigma_at_R(R_perp, R_perp, Sigma1, M[i], concentration, omega_m)
		Deltasigma_1[i] = Deltasigma1

	Size2 = nz * R_perp_bins
	Deltasigma_2 = np.ndarray(Size2)
	Deltasigma_2.shape = (nz,R_perp_bins)

	for i in range(nz):
		xi_mm = ct.xi.xi_mm_at_R(Radii, k_nl, P_k_nl[i])
		Sigma2 = ct.deltasigma.Sigma_at_R(R_perp, Radii, xi_mm, M[0], concentration, omega_m)
		Deltasigma2 = ct.deltasigma.DeltaSigma_at_R(R_perp, R_perp, Sigma2, M[0], concentration, omega_m)
		Deltasigma_2[i] = Deltasigma2

	Size3 = nz * M_bins
	Bias = np.ndarray(Size3)
	Bias.shape = (M_bins,nz)

	#for i in range(M_bins):
	#	for j in range(nz):
	#		Bias[i][j] = ct.bias.bias_at_M(M[i], k, P_k[j], omega_m)


	logM = np.log(M)
	block["deltasigma", "deltasigma_1"] = Deltasigma_1
	block["deltasigma", "deltasigma_2"] = Deltasigma_2
	block["deltasigma", "bias"] = Bias
	block["deltasigma", "m_h"] = M
	block["deltasigma", "lnM"] = logM
	block["deltasigma", "R_perp"] = R_perp
	return 0

def cleanup(config):
	pass
