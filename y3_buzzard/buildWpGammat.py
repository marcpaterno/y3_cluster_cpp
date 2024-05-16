import os
import numpy as np
from cosmosis.datablock import names

from scipy.interpolate import interp1d
from scipy.special import erf
import cluster_toolkit as ct
import massconcen
import bias

from astropy.constants import G
import astropy.cosmology

import hankl
from scipy.fft import ifht

# redshift mean with format [z0, z1, z2, z0, z1, z2, ...]
from setup_bins import zmeans_ij

##################################################################
############# Build Wp(R,z) and Gammat(R,z) ############
# Estimates the clustering of cluster and the shear profiles
# compute a hankel transformation directly from P(k):
#   - Wp_hh/nfw(r)
#   - Shear_hh/nfw(r)
# 
# Returns: 
# WpGammat/
# | Sigma_{hh/nfw}.txt
# | DSigma_{hh/nfw}.txt
# | Bias.tx
# | r.txt
# | z.txt
# | m_h.txt
#
#############################################
# Author: Johnny Esteves
# Created: 20th Dec, 2022
#############################################

cosmo = names.cosmological_parameters
cosmo_fid = astropy.cosmology.FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)

def setup(options):
    section_name = "halo_model"
    #Mpc/h comoving distance, distance on the sky
    R_perp_min = options[section_name,"R_perp_min"]
    R_perp_max = options[section_name,"R_perp_max"]
    R_perp_bins = int(options[section_name,"R_perp_bins"])

    # radii of xi_hm, in Mpc/h
    Radii_min = options[section_name,"Radii_min"]
    Radii_max = options[section_name,"Radii_max"]
    Radii_bins = int(options[section_name,"Radii_bins"])

    #mass (float): Halo mass Msun/h.
    M_min = options[section_name,"M_min"]
    M_max = options[section_name,"M_max"]
    M_bins = int(options[section_name,"M_bins"])

    params_out = R_perp_min, R_perp_max, R_perp_bins, \
    Radii_min, Radii_max, Radii_bins, M_min, M_max, M_bins
    return params_out
    

def execute(block, config):
    section_name = "halo_model"
    R_perp_min, R_perp_max, R_perp_bins, Radii_min, Radii_max,\
    Radii_bins, M_min, M_max, M_bins  = config

    # cosmo parameters
    omega_m = block[cosmo, "omega_m"]
    omega_b = block[cosmo, "omega_b"]
    H0 = block[cosmo, "H0"]
    cosmology = astropy.cosmology.FlatLambdaCDM(H0*100, omega_m, Ob0=omega_b, Tcmb0=2.725)

    # load camb matter-matter power spectrums
    # linear is used for the bias
    k_h = block["matter_power_lin", "k_h"]
    P_k = block["matter_power_lin", "p_k"]
    z_k = block["matter_power_lin", "z"]

    P_k_nl = block["matter_power_nl", "p_k"]
    k_nl = block["matter_power_nl", "k_h"]
    z_nl = block["matter_power_nl", "z"]

    z = z_k

    z_az = block["growth_parameters", "z"]
    daz = block["growth_parameters", "d_z"]
    daz = np.interp(z, z_az, daz)

    # get peak-height values; needs fine-tuning of the k params
    # does not agree with Colossus nor Cluster tool-kit
    # rnu = block["sigma_r", "r"]
    # znu = block["sigma_r", "z"]
    # nu = 1.686/block["sigma_r", "sigma"]

    # compute overdensities; physical not comoving
    rho_c0 = float(cosmology.critical_density0.to('Msun/Mpc^3').value)
    rho_cz = cosmology.critical_density(z).to('Msun/Mpc^3').value
    rho_m = omega_m*rho_c0
    rho_mz = rho_m*(1.+z)**3

    nz = len(z)
    # setup bins
    R_perp = np.logspace(np.log10(R_perp_min), np.log10(R_perp_max), R_perp_bins)
    R_perp = np.append(0, R_perp)
    Radii = np.logspace(np.log10(Radii_min), np.log10(Radii_max), Radii_bins)
    M = np.logspace(np.log10(M_min), np.log10(M_max), M_bins)
    logM = np.log(M)

    # Step 3) Compute Bias (M, Z)
    import time
    t0 = time.time()
    #### Compute peak-height
    Nu0 = ct.peak_height.nu_at_M(M, k_h, P_k[0], omega_m)
    # Involves with 1/da(z)
    # shape (z, mbins)
    nu = np.array([Nu0/daz[i] for i in range(nz)])
    Bias = np.array([bias.bias_at_nu(nui, odelta=200) for nui in nu])
    tf = time.time()-t0
    # print('Bias took: %.3f sec'%tf)

    #### Step 1) Wp
    # compute halo-halo projected correlation function Wp_hh(R)
    # uses hankl fftlog
    # JTA: Wp_hh is the matter-matter correlation function
    # Wp_hh = pk_to_Wp(R_perp, z, k_h, P_k)
    Wp_hh = np.zeros((len(z), R_perp.size),dtype='d')

    # there is no Wp_nf
    # to be removed
    Wp_nfw = Wp_hh #pk_to_Wp(R_perp, z, k_h, pk_nfw)

    #### Step 2) deltaSigma
    ## Cluster Toolkit 2h term computation
    p2h = ct_2hTerm(omega_m, Md=1e14, cd=5, bias=1.0)
    p2h.pk_to_dsigma(R_perp[1:],k_nl,P_k_nl,z_nl)
    Sigma_hh = p2h.Sigma
    dSigma_hh = p2h.dSigma

    # dSigma_hh = np.zeros((nz, Radii.size),dtype='d')

    # compute the NFW projected profile
    # compute NFW Analytical Formula (Wright & Brainerd 2000)
    # using Duffy Mass-Concetration Relation
    # converting M200 to R200 wo the redshift evolution
    # concvec = duffy_concentration_relation(M, z_eff=0.4)
    concvec = child18_mass_concentration(M, 0.4, halo_sample = 'stacked_nfw')

    mm, rr  = np.meshgrid(M, R_perp, indexing='ij')
    dSigma_nfw = deltaSigmaNFW_Analytical(rr, mm, concvec[:,np.newaxis], rho_c=rho_m)
    Sigma_nfw = sigmaNFW_Analytical(rr, mm, concvec[:,np.newaxis], rho_c=rho_m)

    # NOTE: P_k depends on z; mass is a vector
    # Bias has shape (mass.size, z.size)

    # Step 4) Compute the damped halo-halo power spectrum
    # Pk times sqrt(pi)*erf(sigma_photoz * k_nl)/(2*sigma_photoz*k_nl)
    # damped_Pk_nl = compute_damped_pk(sigma_photoz*(1+z), k_nl, P_k_nl)

    # Step 5) shift cosmological scale
    # ratio of the comoving distance dC(cosmo)/dC(cosmo_fid)
    scaleShift, hubbleShift = scaleShiftCosmo(z, cosmology)

    # put into the datablock
    block[section_name, "m_h"] = M
    block[section_name, "lnM"] = logM
    block[section_name, "z"] = z
    block[section_name, "rhoc"] = rho_mz

    # Wp
    block[section_name, "Rp"] = Radii
    block[section_name, "Wp_hh"] = Wp_hh
    block[section_name, "Wp_nfw"] = Wp_nfw

    # deltaSigma and Sigma [Msun/pc^2]
    block[section_name, "r_sigma"] = R_perp
    block[section_name, "DSigma_hh"] = dSigma_hh # already on Msun/pc^2
    block[section_name, "DSigma_nfw"] = dSigma_nfw/1e12
    block[section_name, "Sigma_hh"] = Sigma_hh # already on Msun/pc^2
    block[section_name, "Sigma_nfw"] = Sigma_nfw/1e12
    block[section_name, "concentration"] = concvec
    
    # Bias
    block[section_name, "bias"] = Bias

    # damped Pk
    block[section_name, "scale_shift"] = scaleShift
    block[section_name, "hubble_shift"] = hubbleShift
    block[section_name, "k"] = k_h
    # block[section_name, "Damped_Pk_hh"] = damped_Pk_nl
    
    return 0

def compute_hankel(r, k, pk, mu=0):
    # Hankl means J_l (Bessel 2D) not j_l (Bessel func. 3D)
    si, _res = hankl.P2Wp(k, pk, mu, n=0)
    res = interp1d(si, _res.real)(r)
    return res

def pk_to_Wp(r, z, k, pk):
    """ Compute the Hankel transformation of P(k)

        The projected correlation function Wp(R)

    Args:
        r (array): radial vector
        z (array): redshift vector
        k (array): wave number
        pk (array): power spectrum, shape like (z.size, k.size)

    Returns:
        Wp(R): correlation function
    """
    # start the integration
    Wp = np.zeros((len(z), r.size),dtype='d')
    for i in range(z.size):
        Wp[i] = compute_hankel(r, k, pk[i], mu=0)
    return Wp

## Using Cluster Tool-kit
# The NSIZE of 50 is fast and has an accuracy of 0.1%
# TODO: make sure the dummy variables doesn't affect the result
class ct_2hTerm(object):
    ## following cluster toolkit definition
    RHO_C = 2.77533742639e+11 # Msun/Mpc^3/h^2 (critical density)

    def __init__(self,omega_m=0.3,exclusion=True,NSIZE=50, Md=1e13, cd=4.,bias=1.00):
        self.omega_m = omega_m
        self.NSIZE = NSIZE
        self.Rfix = np.logspace(-3., 3., NSIZE, base=10) #Xi_hm MUST be evaluated to higher than BAO
        self.Md, self.cd = Md, cd # dummy values
        self.bias = bias
        self.exclusion = exclusion
    
    def _pk_to_sigma(self,Rp,k,pk):
        ## TODO: add NFW profile
        xi_mm    = ct.xi.xi_mm_at_r(self.Rfix, k, pk)
        xi_2halo = ct.xi.xi_2halo(self.bias, xi_mm)
        Sigma_mm = ct.deltasigma.Sigma_at_R(Rp, self.Rfix, xi_2halo, self.Md, self.cd, self.omega_m)
        return Sigma_mm
    
    def _to_dsigma(self,Rp,sigma):
        dSigma_mm = ct.deltasigma.DeltaSigma_at_R(Rp, Rp, sigma, self.Md, self.cd, self.omega_m)
        return dSigma_mm

    def pk_to_sigma(self,Rp,k,pk,zvec):
        self.zvec = zvec
        self.Sigma = np.zeros((zvec.size,Rp.size))
        for i in range(zvec.size):
            self.Sigma[i] = self._pk_to_sigma(Rp,k,pk[i])
        pass

    def pk_to_dsigma(self,Rp,k,pk,zvec):
        self.R = Rp
        
        # check if Sigma is computed
        if not hasattr(self,'Sigma'):
            self.pk_to_sigma(Rp,k,pk,zvec)
        
        # convert to delta sigma
        self.dSigma = np.zeros((self.zvec.size,Rp.size))
        for i in range(self.zvec.size):
            sigma = self.Sigma[i]
            if self.exclusion:
                ## ADD NFW
                sigma+= sigmaNFW_Analytical(Rp, self.Md, self.cd, rho_c=self.RHO_C*self.omega_m)/1e12

            self.dSigma[i] = self._to_dsigma(Rp, sigma)
            if self.exclusion:
                ## Subtract NFW
                self.dSigma[i] -= deltaSigmaNFW_Analytical(Rp, self.Md, self.cd, rho_c=self.RHO_C*self.omega_m)/1e12
            
            self.dSigma = np.where(self.dSigma<0.,np.nan,self.dSigma)
        pass

# def _pk_to_dSigma(Rp,k,pk,zvec,omega_m=0.3,NSIZE=50):
#     Md, cd = 1e13, 4. # dummy values
#     Rfix = np.logspace(-3., 3., NSIZE, base=10) #Xi_hm MUST be evaluated to higher than BAO
#     Sigma = np.zeros((zvec.size,Rp.size))
#     dSigma = np.zeros((zvec.size,Rp.size))
#     for i in range(zvec.size):
#         xi_mm    = ct.xi.xi_mm_at_r(Rfix, k, pk[i])
#         xi_2halo = ct.xi.xi_2halo(1.0, xi_mm)
#         Sigma_mm = ct.deltasigma.Sigma_at_R(Rp, Rfix, xi_2halo, Md, cd, omega_m)
#         dSigma_mm = ct.deltasigma.DeltaSigma_at_R(Rp, Rp, Sigma_mm, Md, cd, omega_m)
#         Sigma[i] = Sigma_mm
#         dSigma[i] = dSigma_mm
#     return Sigma, dSigma

# JTA: There should be a reference to the literature on this equation
# def pk_to_dSigma(r, z, k, pk):
#     """ Compute the second-order Hankel transformation of P(k)

#         The projected density profile dSigma(R)  (aka deltaSigma(R) )

#     Args:
#         r (array): radial vector
#         z (array): redshift vector
#         k (array): wave number
#         pk (array): power spectrum, shape like (z.size, k.size)

#     Returns:
#         dSigma(R): projected density profile
#     """
#     # start the integration
#     dS = np.zeros((len(z), r.size),dtype='d')
#     for i in range(z.size):
#         dS[i] = compute_hankel(r, k, pk[i], mu=2)
#     return dS

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
    # double loop is faster
    # I know it is weird
    for i in range(mass.size):
        for j in range(zsize):
            try:
                Bias[i][j] = ct.bias.bias_at_M(mass[i], k, P_k[j], omega_m)
            except:
                print('Bias error')
                Bias[i][j] = 0.
    return Bias

def compute_concentration(z, mass, mstar, z_mstar, 
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
        conc[i] = massconcen.c_from_m200(mass[i], z, omega_m, omega_b, 
                                          sigma8, h0, mstar, z_mstar)
    return conc

def scaleShiftCosmo(znew, cosmo, eps=1e-9):
    """Scale Shift Cosmology

    To adapt to a new cosmology we can re-scale the distances by taking
    into account the fiducial cosmolgoy.
    
    scaleShiftCosmo is basically the ratio of the comoving distance
    and the hubble factor
        scaleShiftCosmo = (dist_c/h) / (dist_c_fid/h_fid)
    
    the fiducial cosmology was set to Omega_m = 0.3 and H0 = 70

    Args:
        znew (array): redshift vector
        cosmo (astropy.cosmology): astropy cosmology

    Returns:
        scale_shift: scale shift factor = ratio of comoving distances/H(z)
        hubble_shift: the hubble flow is the shift on the parallel direction of the los
    """
    # cosmological quantities
    # h0 = block[cosmo, "h0"]

    # fiducical cosmology
    Hz_fid = cosmo_fid.H(znew).value
    dc_fid = cosmo_fid.comoving_distance(znew).value

    # current cosmology
    Hz = cosmo.H(znew).value
    dc = cosmo.comoving_distance(znew).value

    # scale shift
    scale_shift = dc/(dc_fid+eps)
    scale_shift[0] = 1.

    # hubble shift
    hubble_shift = Hz/Hz_fid

    return scale_shift, hubble_shift

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

from scipy.special import sici
# def pk_nfw_profile(k_h, m200, rho_c=1., omega_m=0.3, z_eff=0.4):
#     c = duffy_concentration_relation(m200, z_eff=z_eff)
#     rho_cz = rho_c*(omega_m*(1+z_eff)**3+(1-omega_m))
#     r_virial = convert_m200_to_r200(m200, rho_cz)
#     rho_k = (m200/(rho_c*omega_m))*normalized_halo_profile(k_h, r_virial, c)
#     return rho_k

def duffy_concentration_relation(m_h, z_eff=0.4):
    a_eff = 1/(1+z_eff)
    m_h_pivot = 2e12
    return 7.85*np.power(m_h/m_h_pivot,-0.081)*np.power(a_eff,0.71)

def convert_m200_to_r200(m200,rho,odelta=200):
    rv = np.power(m200*3.0/(4.0*np.pi*rho*odelta),1.0/3.0)
    return rv

def normalized_halo_profile(k_h, r_virial, c):
    """Compute the normalized halo profile function in Fourier space; :math:`u(k|m)`

    For details of the available profile parametrizations, see the class description.

    Note that the function returns unity for :math:`k \\leq 0`.

    Args:
        k_h (np.ndarray): Wavenumber in h/Mpc units.
        r_virial (np.ndarray): Virial radius in Mpc/h units.
        c (np.ndarray): Halo concentration parameter; :math:`c = r_\mathrm{virial}/r_\mathrm{scale}`.

    Returns:
        np.ndarray: Normalized halo profile :math:`u(k|m)`

    copied from effective halos package
    """
    r_scale = r_virial/c # in Mpc/h units
    ks = k_h*r_scale

    sici1 = sici(ks)
    sici2 = sici(ks*(1.+c))
    f1 = np.sin(ks)*(sici2[0]-sici1[0])
    f2 = np.cos(ks)*(sici2[1]-sici1[1])
    f3 = np.sin(c*ks)/(ks*(1.+c))
    fc = np.log(1.+c)-c/(1.+c)
    return (f1+f2-f3)/fc

### Compute \Sigma Analyticaly
def f_greater_than(x):
    res = 1. / (x**2 - 1.0)
    res *=(1- 2.0 / np.sqrt(x**2 - 1.0) * np.arctan(np.sqrt((x - 1.0) / (1.0 + x))))
    return res

def f_less_than(x):
    res = 1. / (x**2 - 1.0)
    res *=(1- 2.0 / np.sqrt(1.0 - x**2) * np.arctanh(np.sqrt((1.0 - x) / (1.0 + x))))
    return res
        
def f_nfw(x, eps=1e-9):
    """
    Analytical normalized Sigma profile
    Wright & Brained 2000

    Args:
        x (array): Rp/Rs where Rs is the scale radius, Rs = R200/c
    """
    if (isinstance(x,float))or(isinstance(x,int)): 
        x = np.array([x])
        
    res = 1/3.*np.ones_like(x)

    ix = np.where(x <= 1-eps)[0]
    res[ix] = f_less_than(x[ix])
    
    ix = np.where(x>=1+eps)[0]
    res[ix] = f_greater_than(x[ix])
    return res

def sigmaNFW_Analytical(radii, m200, c, rho_c=0.0):
    """Analytical NFW Surface Mass Density  (Eqn 14, Wright & Brained 2000)

    Args:
        radii (array): projected radii
        m200 (array): M200 critical mass in solar masses unit
        z_eff (float): effective redshift to concentration and critical density values
    """
    r_virial = convert_m200_to_r200(m200, rho_c)
    rs = r_virial/c

    # eqn 3 
    deltac = (200./3.) * c**3/(np.log(1+c)-c/(1+c))

    # eqn 15 an 16
    # interpolation avoid numerical issues
    xvec = np.logspace(-3.,4.,1000)
    xvec = np.append(0, xvec)
    fx = interp1d(xvec, f_nfw(xvec))(radii/rs)
    return 2*rs*deltac*rho_c*fx

### Compute \Delta \Sigma Analyticaly
def g_less_than(x,core=4e-3):
    # below core the solution fails numereically
    x = np.where(x<core, core, x)
    term1 = 8.0*np.arctanh(np.sqrt((1.0-x)/(1.0+x)))/(x**2*np.sqrt(1.0-x**2))
    term2 = 4.0/x**2 * np.log(x/2.0)
    term3 = -2.0/(x**2-1.0)
    term4 = 4.0*np.arctanh(np.sqrt((1.0-x)/(1.0+x)))/((x**2-1.0)*np.sqrt(1.0-x**2))
    return term1 + term2 + term3 + term4

def g_greater_than(x, xmax=1e9):
    # above xmax the function achieves machine precision
    x = np.where(x>xmax, xmax, x)
    term1 = 8.0*np.arctan(np.sqrt((x-1.0)/(1.0+x)))/(x**2*np.sqrt(x**2-1.0))
    term2 = 4.0/x**2 * np.log(x/2.0)
    term3 = -2.0/(x**2-1.0)
    term4 = 4.0*np.arctan(np.sqrt((x-1.0)/(1.0+x)))/((x**2-1.0)**(3.0/2.0))
    return term1 + term2 + term3 + term4

def g_nfw(x, eps=1e-9):
    """gNFW Eqn 15 and 16 Wright & Brained 2000

    Analytical normalized shear/deltaSigma profile

    Args:
        x (array): Rp/Rs where Rs is the scale radius, Rs = R200/c
    """
    res = np.zeros_like(x)
    ix = np.where(np.abs(x-1) <= eps)[0]
    res[ix] = 10./3. + 4*np.log(1/2.) 
    
    ix = np.where(x <= 1-eps)[0]
    res[ix] = g_less_than(x[ix])
    
    ix = np.where(x>=1+eps)[0]
    res[ix] = g_greater_than(x[ix])
    return res

def convert_m200_to_r200(m200,rho,odelta=200):
    rv = np.power(3*m200/(4.0*np.pi*rho*odelta),1.0/3.0)
    return rv

def deltaSigmaNFW_Analytical(radii, m200, c, rho_c=0.0):
    """Analytical NFW Surface Mass Density  (Eqn 14, Wright & Brained 2000)

    Args:
        radii (array): projected radii
        m200 (array): M200 critical mass in solar masses unit
        z_eff (float): effective redshift to concentration and critical density values
    """
    r_virial = convert_m200_to_r200(m200, rho_c)
    rs = r_virial/c

    # eqn 3 
    deltac = (200./3.) * c**3/(np.log(1+c)-c/(1+c))

    # eqn 15 an 16
    # interpolation avoid numerical issues
    xvec = np.logspace(-3.,4.,1000)
    xvec = np.append(0, xvec)
    gx = interp1d(xvec, g_nfw(xvec))(radii/rs)
    return rs*deltac*rho_c*gx

def compute_bias_camb(rnu, znu, nu, massvec, zvec, rhomz):
    # loop over redshift
    # coult it be improved; but it is fast enough
    
    Bias = np.zeros((massvec.size, zvec.size))
    for i,z in enumerate(zvec):
        # slice nu(radii,z)
        iz = np.interp(z,znu,np.arange(znu.size)).astype(int)
        nui = nu[iz]

        # convert mass to radii_lagrangian
        # Child et al. 2018 Eqn 13
        rhom = rhomz[i]
        rvec = convert_m200_to_r200(massvec, rhom, odelta=1)
        
        # interpolate with the new radii vector
        _nui = np.exp(interp1d(np.log(rnu), np.log(nui))(np.log(rvec)))
        Bias[:,i] = bias.bias_at_nu(_nui, odelta=200)
    return Bias

def child18_mass_concentration(M200c, z, halo_sample = 'stacked_nfw'):
	"""
	The concentration mass relation model of Child et al 2018.

	Parameters
	-----------------------------------------------------------------------------------------------
	M200c: array_like
		Halo mass in :math:`M_{\odot}/h`; can be a number or a numpy array.
	z: float
		Redshift
	halo_sample: str
		Can be ``individual_all`` (default), ``individual_relaxed`` (the mean concentration of 
		individual, relaxed halos), ``stacked_nfw`` (the stacked profile with with an NFW profile), 
		and ``stacked_einasto`` (the stacked profile with with an Einasto profile).
		
	Returns
	-----------------------------------------------------------------------------------------------
	c: array_like
		Halo concentration; has the same dimensions as ``M``.
	"""

	if halo_sample == 'individual_all':
		m = -0.10
		A = 3.44
		b = 430.49
		c0 = 3.19
	elif halo_sample == 'individual_relaxed':
		m = -0.09
		A = 2.88
		b = 1644.53
		c0 = 3.54
	elif halo_sample == 'stacked_nfw':
		m = -0.07
		A = 4.61
		b = 638.65
		c0 = 3.59
	elif halo_sample == 'stacked_einasto':
		m = -0.01
		A = 63.2
		b = 431.48
		c0 = 3.36
	else:
		raise Exception('Unknown halo sample for child18 concentration model, %s.' % (halo_sample))

	Mstar = peakHeight_nonLinearMass(z)
	M_MT = M200c / (Mstar * b)	
	c200c = c0 + A * (M_MT**m * (1.0 + M_MT)**-m - 1.0)
	return c200c

def peakHeight_nonLinearMass(z):
    """peakHeight_nonLinearMass 
    The non-linear mass. E.g. eqaution 13 in Child et al. 2018
    Mstar is log-linear with redshift. 
    log(Mstar) = 12.5 - 1.5*z

    You can check with Colossus
    How to generate the data:
    from colossus.cosmology import cosmology
    from colossus.lss.peaks import nonLinearMass
    cosmology.setCosmology('WMAP7-only')
    z = np.linspace(0, 2, 100)
    Mstar = nonLinearMass(z)
    """
    return 10**(12.5-1.5*z)


def cleanup(config):
    pass

