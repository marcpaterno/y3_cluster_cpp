import os
import numpy as np
from cosmosis.datablock import option_section, names

from scipy.interpolate import interp1d
from scipy.special import erf
import cluster_toolkit as ct
import massconcen

import hankl

# redshift mean with format [z0, z1, z2, z0, z1, z2, ...]
from setup_bins import zmeans_ij

##################################################################
############# Build Xi(r,z), Sigma(r,z) and bias(M,z) ############
# Estimates the hala-halo and nfw correlation functions (CF),
# projected halo density (Sigma) and the bias:
#   - Xi_hh/nfw(r)
#   - Sigma_hh/nfw(r)
#   - b(M)
#   - P_hh_damp(k)
# 
# Estimates Sigma(r) by the sky projection of the CFs:
#     \Sigma(r) = \rho_m \int Xi(\sqrt{R^2+z^2}) dz
#
# see https://cluster-toolkit.readthedocs.io/en/latest/source/projected_density_profiles.html
#
# We also compute the bias of the matter-matter CF:
#     \Xi_hh = b(M) \Xi_mm
#
# see https://cluster-toolkit.readthedocs.io/en/latest/source/halo_bias.html
#
# As a helper we compute the damped matter-matter power spectrum:
#      P_hh_damp(k_s, k_z) \propto efc(\sigma_z * k) P_hh(k)/\sigma_z * k
# 
# the \sigma_{z,0} is the cluster photo-z redshift error. 
# See equation (8) https://arxiv.org/abs/0801.3485.
#
# The halo-halo CF is the hankel transformation
# of the halo-halo power spectrum. To perform
# this transformation we use FFTLog described
# in hankl ...
# 
# The Xi_nfw(r) is basically:
#  \xi_nfw(r) = (\rho_nfw(r)}/\Omega_m\rho_crit) - 1
# 
# Returns: 
# corrFunctions/
# | Xi_hh/nfw.txt
# | Pk_hh_damp.txt
# | Sigma_hh/nfw.txt
# | r_xi.txt
# | r_sigma.txt
# | z.txt
# | m_h.txt
#
#############################################
# Author: Johnny Esteves
# Created: 17th Nov, 2022
#############################################

cosmo = names.cosmological_parameters


def setup(options):
    section = option_section
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

    # photo-z cluster error, used in the damping factor
    sigma_photoz= options[section,"sigma_photoz"]

    params_out = R_perp_min, R_perp_max, R_perp_bins, \
    Radii_min, Radii_max, Radii_bins, M_min, M_max, M_bins, sigma_photoz
    return params_out
    

def execute(block, config):
    R_perp_min, R_perp_max, R_perp_bins, Radii_min, Radii_max,\
    Radii_bins, M_min, M_max, M_bins, sigma_photoz  = config

    # cosmo parameters
    omega_m = block[cosmo, "omega_m"]
    omega_b = block[cosmo, "omega_b"]
    sigma8 = block[cosmo, "sigma8_input"]
    h0 = block[cosmo, "h0"]

    # load camb matter-matter power spectrums
    # linear is used for the bias
    k = block["matter_power_lin", "k_h"]
    P_k = block["matter_power_lin", "p_k"]
    z_k = block["matter_power_lin", "z"]

    k_nl = block["matter_power_nl", "k_h"]
    P_k_nl = block["matter_power_nl", "p_k"]
    z_k_nl = block["matter_power_nl", "z"]
    
    z = block["matter_power_nl", "z"]
    nz = len(z)

    # used in the mass-concetration relation
    # load charcteristic mass
    mstar = block["mstar", "mstar"]
    mstar_z = block["mstar", "z"]

    # setup bins
    R_perp = np.logspace(np.log10(R_perp_min), np.log10(R_perp_max), R_perp_bins)
    Radii = np.logspace(np.log10(Radii_min), np.log10(Radii_max), Radii_bins)
    M = np.logspace(np.log10(M_min), np.log10(M_max), M_bins)
    logM = np.log(M)

    # compute NFW concentration
    # assume mass-concentration relation
    zp = 0.4 # assuming at redshift 0.4s
    cM = compute_conentration(zp, M, mstar, mstar_z,
                              CosmoParams=(omega_m, omega_b, sigma8, h0))

    #### Step 1) Correlation Functions

    # compute halo-halo correlation function Xi_hh(r)
    # uses cluster-toolkit
    # TODO switch to FFTLog
    Xi_hh, s = pk_to_Xi(Radii, z, k_nl, P_k_nl)

    # the mass concentration relation is a function of redshift, but the way Y. Zhang designed 
    # this, the 1-halo term is independent of redshift, which is computationally efficient.
    # We'll have to assume a mean redshift and compute the c from M at that z.

    # compute nfw correlation function Xi_nfw(r)
    Xi_nfw = compute_Xi_nfw(Radii, M, cM, omega_m)

    # Do 2d projection
    Wp_hh = xi_to_wp(s, Xi_hh, pimax=100)
    Wp_nfw = xi_to_wp(Radii, Xi_nfw, pimax=100)

    #### Step 2) Surface Mass Density Profiles - Sigma [h Msun/pc^2]

    # compute the halo-halo projected sigma profile
    Sigma_hh = compute_Sigma_hh(R_perp, Radii, Xi_nfw, omega_m)

    # compute the NFW projected profile
    Sigma_nfw = compute_Sigma_nfw(R_perp, M, cM, omega_m)

    # Step 3) Compute Bias (M, Z)
    Bias = compute_bias(M, k, P_k, omega_m)
    # NOTE: P_k depends on z; mass is a vector
    # Bias has shape (mass.size, z.size)

    # Step 4) Compute the damped halo-halo power spectrum
    # Pk times sqrt(pi)*erf(sigma_photoz * k_nl)/(2*sigma_photoz*k_nl)
    damped_Pk_nl = compute_damped_pk(sigma_photoz*(1+z), k_nl, P_k_nl)

    # Step 5) shift cosmological scale
    # ratio of the comoving distance dC(cosmo)/dC(cosmo_fid)
    scaleShift, hubbleShift = scaleShiftCosmo(z, block, h0=h0, omega_m=omega_m)

    # put into the datablock
    block["correlationFunction", "m_h"] = M
    block["correlationFunction", "lnM"] = logM
    block["correlationFunction", "z"] = z

    # Xi
    block["correlationFunction", "r_xi"] = s
    block["correlationFunction", "Xi_hh"] = Xi_hh
    block["correlationFunction", "Xi_nfw"] = Xi_nfw

    # Wp
    block["correlationFunction", "Rp"] = Radii
    block["correlationFunction", "Wp_hh"] = Wp_hh
    block["correlationFunction", "Wp_nfw"] = Wp_nfw

    # Sigma
    block["correlationFunction", "r_sigma"] = R_perp
    block["correlationFunction", "Sigma_nfw"] = Sigma_nfw
    block["correlationFunction", "Sigma_hh"] = Sigma_hh
    
    # Bias
    block["correlationFunction", "bias"] = Bias

    # damped Pk
    block["correlationFunction", "scale_shift"] = scaleShift
    block["correlationFunction", "hubble_shift"] = hubbleShift
    block["correlationFunction", "k"] = k
    block["correlationFunction", "Damped_Pk_hh"] = damped_Pk_nl

    # Debug
    pk_nl_in = pk_nl[25]
    s, xi_nl = hankl.P2xi(k_nl, pk_nl_in, l=0)
    xi_nl_ct = ct.xi.xi_mm_at_r(s, k_nl, pk_nl_in)
    np.savez('./pk_nl.npz',k=k_nl,pk=pk_nl_in)
    np.savez('./xi_hankl_nl.npz',r=s,pk=xi_nl)
    np.savez('./xi_ct_nl.npz',r=s,pk=xi_nl_ct)
    
    return 0
        
def pk_to_Xi(r, z, k, Pk):
    """ Compute the hankel transformation of P(k)

        The correlation function Xi(r)

    Args:
        r (array): radial vector
        z (array): redshift vector
        k (array): wave number
        Pk (array): power spectrum, shape like (z.size, k.size)

    Returns:
        Xi(r): correlation function
    """
    # start the integration
    Xi = np.zeros((len(z), k.size))
    for i in range(z.size):
        # using cluster toolkit
        # switch to FFTlog
        # Xi[i] = ct.xi.xi_mm_at_r(r, k, Pk[i])
        
        si, xii = hankl.P2xi(k, Pk[i], l=0)
        # Xi[i] = interp1d(si, xii.real)(r)
        Xi[i] = xii

    return Xi, si

def compute_Xi_nfw(r, mass, c, omega_m=0.3):
    """ NFW Profile Correlation Function

    Computes Xi_nfw using cluster-tool kit implementation

    Args:
        r (array): radius
        mass (array): halo 200 mass
        c (array): concentration, same size of halo mass
        omega_m (float, optional): matter fraction. Defaults to 0.3.

    Returns:
        Xi (ndarray): Correlation function with shape (mass.size, r.size)
    """ 
    Xi = np.zeros((mass.size, r.size))
    for i in range(mass.size):
        Xi[i] = ct.xi.xi_nfw_at_r(r, mass[i], c[i], omega_m)
    return Xi

def compute_Sigma_hh(R, r, Xi, omega_m=0.3,
                     dummy_conc=4., dummy_mass=3.2e14):
    """ Surface Mass Density Profile (Sigma)

    For any correlation function computes the surf. mass
    density profile by computing the projection on the sky.

    Computes Sigma using cluster-tool kit implementation

    Note: the dummy variables are just there to keep the 
    args consistent in cluster-tool kit. they can be any
    arbritray number.
    https://cluster-toolkit.readthedocs.io/en/latest/source/projected_density_profiles.html

    Args:
        R (array): projected radius
        r (array): 3d radius
        Xi (ndarray): correlation function (r)
        omega_m (float, optional): matter fraction. Defaults to 0.3.

    Returns:
        Sigma (ndarray): Projected NFW Profile, units h Msun/pc^2,
        (mass.size, r.size)
    """ 
    zsize = Xi.shape[0]
    Sigma = np.zeros((zsize, R.size))

    for i in range(zsize):
        Sigma[i] = ct.deltasigma.Sigma_at_R(R, r, Xi[i],
                                           dummy_mass, dummy_conc, omega_m)
    return Sigma

def compute_Sigma_nfw(R, mass, c, omega_m=0.3):
    """ Projected NFW Profile (Sigma_nfw)

    Computes Sigma_nfw using cluster-tool kit implementation

    Args:
        R (array): projected radius
        mass (array): halo 200 mass
        c (array): concentration, same size of halo mass
        omega_m (float, optional): matter fraction. Defaults to 0.3.

    Returns:
        Sigma (ndarray): Projected NFW Profile, units h Msun/pc^2,
        (mass.size, r.size)
    """ 
    Sigma = np.zeros((mass.size, R.size))
    for i in range(mass.size):
        Sigma[i] = ct.deltasigma.Sigma_nfw_at_R(R, mass[i], c[i], omega_m)
    return Sigma

def compute_damped_pk(sigma_photoz, k, Pk):
    """ Power Spectrum Times a Photo-z Damping Factor

    Accounts for the cluster photo-z error on Pk

    Args:
        sigma_photoz (array or float): photo-z cluster error
        k (array): wave number
        Pk (ndarray): power spectrum

    Returns:
        damped_pk: shape like pk
    """
    if len(sigma_photoz)==0:
        sigma_photoz = sigma_photoz*np.ones_like(Pk[:,0])
    
    const = np.sqrt(np.pi)/2.
    # multiply the power spectrum by the damping factor
    damped_Pk = np.zeros_like(Pk)
    for i in range(Pk.shape[0]):
        dampingFactor = const*erf(sigma_photoz[i] * k)/(sigma_photoz[i] * k)
        damped_Pk[i]  = dampingFactor * Pk[i]

    return damped_Pk

def compute_bias(mass, k, P_k, omega_m=0.3):
    """ Compute Matter-Matter Pk Bias

        Xi_hh = bias Xi_mm

    Args:
        mass (array): halo mass
        k (array): wave number
        P_k (ndarray): linear power spectrum as a function of redshift
        omega_m (float, optional): matter fraction. Defaults to 0.3.
    
    Returns:
        bias (ndarray): (mass.size, z.size)

    """
    zsize = P_k.shape[0]
    Bias = np.zeros((mass.size, zsize))
    for j in range(zsize):
            Bias[:, j] = ct.bias.bias_at_M(mass, k, P_k[j], omega_m)
    return Bias

def compute_conentration(z, mass, mstar, z_mstar, 
                         CosmoParams=(0.3, 0.05, 0.82, 0.7)):
    """Concetration-Mass Relation 

    Args:
        z (float): halo redshift 
        mass (array): critical halo mass (200)
        mstar (array): ...tbd - ask Jim
        z_mstar (array): interpolation vector of mstar
        CosmoParams (tuple, optional): cosmological values:
         omega_m, omega_b, sigma8, h0. Defaults to (0.3, 0.05, 0.82, 0.7).
    """
    omega_m, omega_b, sigma8, h0 = CosmoParams
    conc = np.zeros(mass.size)
    for i in range(mass.size):
        conc[i] = massconcen.c_from_m200 (mass[i], z, omega_m, omega_b, 
                                          sigma8, h0, mstar, z_mstar)
    return conc

from astropy.cosmology import FlatLambdaCDM
cosmo_fid = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)

def scaleShiftCosmo(znew, block, h0=0.7, omega_m=0.3):
    """Scale Shift Cosmology

    To adapt to a new cosmology we can re-scale the distances by taking
    into account the fiducial cosmolgoy.
    
    scaleShiftCosmo is basically the ratio of the comoving distance
    and the hubble factor
        scaleShiftCosmo = (dist_c/h) / (dist_c_fid/h_fid)
    
    the fiducial cosmology was set to Omega_m = 0.3 and H0 = 70

    Args:
        znew (array): redshift vector
        block (dict): data block
        h0 (float, optional): small hubble constant. Defaults to 0.7.

    Returns:
        scale_shift: scale shift factor
    """
    # cosmological quantities
    # h0 = block[cosmo, "h0"]

    z_dc = block["distances",'z']
    dc = block["distances",'d_a']*(1+z_dc) # comoving distance Mpc

    z = block["growth_parameters","z"]
    hz = block["growth_parameters","h"]

    # fiducical cosmology
    H0_fid = cosmo_fid.H(z_dc).value
    dc_fid = cosmo_fid.comoving_distance(z_dc).value

    # scale shift
    scale_shift_vec = dc/dc_fid
    scale_shift_vec[0] = 1.

    # the hubble flow is th shift on the parallel direction of the los
    H0_vec = np.interp(z_dc, z, hz)
    hubble_shift_vec = H0_fid/H0_vec

    # interpolate for the new redshift
    scale_shift_perp = interp1d(z_dc[1:], scale_shift_vec[1:], fill_value='extrapolate')(znew)
    # avoids dc = 0

    hubble_shift = np.interp(znew, z_dc, hubble_shift_vec)
    return scale_shift_perp, hubble_shift

def xi_to_wp(r, xi, pimax=100):
    """ Convert Xi(r) to Wp(R)

        Project Xi(r) onto the sky plane by perfoming an Abell Transform.

    Args:
        r (array): radial vector
        z (array): redshift vector
        pimax (float): max values of the parallel direction
        Wp (array): shape like (z or M size, r size)

    Returns:
        Wp(R): projected correlation function
    """
    # start the integration
    Wp = np.zeros_like(xi)
    for i, xii in enumerate(xi):
        Wp[i] = compute_abell_transform(r, xii, pimax=pimax, int_func=np.trapz)
    return Wp    


def compute_abell_transform(r, fr, pimax=100., int_func=np.trapz):
    # Computes Wp(R) = \int \xi(\sqrt{R^2+\pi^2}) d\pi
    #
    # The equation above can be re-written to:
    #
    # Wp(R) = \int (dr/\sqrt{r^2-R^2}) r \xi(r).
    #
    # Calculation of the integral  used in Abel transform
    #          r_max
    #         ⌠
    #         ⎮     r \xi(r)
    #         ⎮ ────────────── dr
    #         ⎮    ___________
    #         ⎮   ╱  2   2
    #         ⎮ ╲╱  r - R
    #         ⌡
    #         R
    # Where r_max=\sqrt{R^2+\pi_{max}^2}
    # Code based on: PyAbel package, direct method.
    # https://github.com/PyAbel/PyAbel/blob/master/abel/direct.py
    # TODO: Implement the C++ version.
    
    # Switch to the r coordinates.
    fr *= 2*r

    if is_uniform_sampling(r):
        int_opts = {'dx': abs(r[1] - r[0])}
    else:
        int_opts = {'x': r}

    # Initializiate output
    out = np.zeros(fr.size)

    # Initialize grid
    R, rr = np.meshgrid(r, r, indexing='ij')
    i_vect = np.arange(len(r), dtype=int)
    II, JJ = np.meshgrid(i_vect, i_vect, indexing='ij')

    # Set integration limits
    # only for R < r <= rmax 
    mask = (II < JJ)
    mask &= (rr <= np.sqrt(R**2 + pimax**2))

    I_sqrt = np.zeros(R.shape)
    I_sqrt[mask] = np.sqrt((rr**2 - R**2)[mask])

    I_isqrt = np.zeros(R.shape)
    I_isqrt[mask] = 1./I_sqrt[mask]

    # create a mask that just shows the first two points of the integral
    mask2 = ((II > JJ-2) & (II < JJ+1))
    
    # perfom the integration
    Integrand = fr[None, :] * I_isqrt  # set up the integral
    out = int_func(Integrand, axis=1, **int_opts)  # take the integral
    # correct for the extra triangle at the start of the integral
    out = out - 0.5*int_func(Integrand*mask2, axis=1, **int_opts)

    return out

def is_uniform_sampling(r):
    """
    Returns True if the array is uniformly spaced to within 1e-13.
    Otherwise False.
    """
    dr = np.diff(r)
    ddr = np.diff(dr)
    return np.allclose(ddr, 0, atol=1e-13)

def cleanup(config):
    pass

