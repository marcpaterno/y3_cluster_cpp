Simplified GPU Modules
=================================================

These modules implement cluster lensing observables using direct MOR integration.

MODULES
-------

1. numberCountsSimple.cu
   - Computes number counts: N = int n(M,z) * P(lo|M,z) * K(zo|z)
   - Grid: (zo_low, zo_high) pairs for redshift bins
   - Output: Total counts in each bin

2. sigma1hSimple.cu
   - Computes N * <Sigma_1h(R)>: 1-halo surface density weighted by counts
   - Grid: (zo_low, zo_high, radii)
   - Divide by numberCountsSimple to get bin-averaged <Sigma_1h>

3. sigmaTotSimple.cu
   - Computes N * <Sigma_tot(R)>: Total (1h + b*2h) surface density
   - Grid: (zo_low, zo_high, radii)
   - Uses sig_max.cuh which computes Sigma_1h + b(M,z) * Sigma_2h

4. shearTotSimple.cu
   - Computes N * <gamma_t(R)>: Total tangential shear
   - Grid: (zo_low, zo_high, radii)
   - Uses gamma_max.cuh which computes DSigma_tot * Sigma_crit_inv

5. numberCountsFull.cu  [NEW - Full 4D Integration]
   - Computes number counts with FULL physics:
     N = int dlob int dltr int dz int dlnM
         Omega(z) * dV/dOmega/dz * n(M,z) * P_HOD(ltr|M,z) * P_EMG(lob|ltr,z) * K(zo|z)
   - 4D integration over (lambda_ob, lambda_tr, z, lnM)
   - Includes EMG observation scatter P(lob|ltr,z) from emg_des_t.cuh
   - Grid: (lob_bin_low, lob_bin_high, zo_low, zo_high) for richness+redshift bins
   - Purpose: Validate GPU pipeline against CPU sel_function approach
   - Note: 4D integration is slower but matches full physics

HELPER MODULES
--------------

plob_ltr_loader.py
   - Loads EMG coefficients from data/prj_params/plob_ltr_params.npz
   - Required by numberCountsFull.cu
   - Writes to datablock section 'plob_ltr_params'

DIFFERENCES FROM sel_function PIPELINE
--------------------------------------

The standard pipeline (NumCountsSel, Sigma1hSel, etc.) uses:
  - sel_function module to pre-tabulate S_ij(lnM, z)
  - S_ij encodes both MOR P(lo|M,z) and richness kernel K_i(lo|lt,z)
  - Integration is over (zt, lnM) only, S_ij absorbs the richness dimension

These simplified modules:
  - Integrate directly over (lo, zt, lnM)
  - Use MOR P(lo|M,z) from mor_des_log_t.cuh
  - Use photo-z kernel K(zo|z) from int_zo_zt_des_t.cuh
  - Do NOT include richness selection effects (projection bias, etc.)
  - Suitable for validation against simulations or simplified models

REQUIRED DATABLOCK SECTIONS
---------------------------

For simple modules (numberCountsSimple, sigma1hSimple, etc.):
- cosmological_parameters: omega_M, omega_nu
- cluster_abundance: MOR parameters, hmf_s, hmf_q
- mass_function: m_h, z, dndlnmh
- photoz: sigma_z, z
- haloModel: r_sigma, lnM, z, Sigma_nfw, Sigma_hh, DSigma_nfw, DSigma_hh,
             bias, sigmaCritInv

For numberCountsFull (4D integration):
- cosmological_parameters: omega_M, omega_nu
- cluster_abundance: MOR parameters (mor_logMmin, mor_logRatio, mor_alpha,
                     mor_sigma, mor_epsilon, z_mor_pivot), hmf_s, hmf_q
- mass_function: m_h, z, dndlnmh
- plob_ltr_params: z, a_mu, b_mu, a_sig, b_sig, a_tau, b_tau, a_fprj, b_fprj
                   (loaded by plob_ltr_loader.py)

EXAMPLE INI CONFIGURATION
-------------------------

[numberCountsSimple]
file = ${Y3_CLUSTER_CPP_DIR}/release-build/src/modules/simple_gpu/numberCountsSimple.so

; Integration bounds
lo_low = 10.0
lo_high = 200.0
zt_low = 0.1
zt_high = 0.9
lnm_low = 31.0
lnm_high = 36.0

; Grid points (wall-of-numbers format)
zo_low = 0.2 0.35 0.5
zo_high = 0.35 0.5 0.65

; Integration settings
eps_rel = 1e-3
eps_abs = 1e-10
max_eval = 1.0e6
algorithm = pagani


EXAMPLE INI FOR numberCountsFull (4D)
-------------------------------------

[plob_ltr_loader]
file = ${Y3_CLUSTER_CPP_DIR}/src/modules/simple_gpu/plob_ltr_loader.py

[numberCountsFull]
file = ${Y3_CLUSTER_CPP_DIR}/release-build/src/modules/simple_gpu/numberCountsFull.so

; 4D integration bounds: (lob, ltr, zt, lnM)
lob_int_low = 5.0
lob_int_high = 250.0
ltr_low = 1.0
ltr_high = 300.0
zt_low = 0.1
zt_high = 0.8
lnm_low = 31.0
lnm_high = 36.0

; Grid points: richness and redshift bins
; DES Y3 bins: lambda=[20,30,45,60,200], z=[0.20,0.35,0.50,0.65]
lob_bin_low = 20.0 20.0 20.0 30.0 30.0 30.0
lob_bin_high = 30.0 30.0 30.0 45.0 45.0 45.0
zo_low = 0.20 0.35 0.50 0.20 0.35 0.50
zo_high = 0.35 0.50 0.65 0.35 0.50 0.65

; Integration settings (4D needs more evaluations)
eps_rel = 1e-2
eps_abs = 1e-10
max_eval = 5.0e6
algorithm = pagani
use_cartesian_product = false


NEW CUDA MODELS
---------------

emg_des_t.cuh (src/models/)
   - EMG (Exponentially Modified Gaussian) kernel for P(lambda_ob | lambda_tr, z)
   - Implements the Costanzi observation scatter model
   - Parameters (mu, sigma, tau, fprj) are z-dependent via spline interpolation
   - Methods:
     * get_params(ltr, z, mu, sigma, tau, fprj) - compute EMG parameters
     * cdf(lob, ltr, z) - CDF of the EMG distribution
     * operator()(lob, ltr, z) - PDF of the EMG distribution
