import numpy as np
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.interpolate import interp1d

def setup(options):
    section = option_section
    # data vector 
    profiles = np.genfromtxt('data_files/wl_data_vector.txt')
    covmat = np.genfromtxt('data_files/wl_cov.txt')
    ncs = np.genfromtxt('data_files/nc.txt', delimiter=',')
    nc_covmat = np.genfromtxt('data_files/Cov_ij_bestfit_DESY1_105.txt')
    ncs = ncs[0, :]
    nc_covmat = nc_covmat[:4, :4]
    
    # model module     
    '''
    profile_blockname = options[section,"profile_blockname"]
    NCs_blockname = options[section,"NCs_blockname"]
    snapshot_zs = options[section, "snapshot_zs"]
    snapshot_lnm = options[section, "lnm_low"]
    snapshot_radii = options[section, "radii"]
    '''

    rad = options[section, "radii"]
    Redges = np.logspace(np.log10(0.0323), np.log10(30.), num=15+1)
    return profiles, covmat, ncs, nc_covmat, rad, Redges

def execute(block, config):

    #read in the model values at this sample point. 
    profiles_miscent_model = block['SigmaMiscentY1CUDAIntegrand', "vals"]
    profiles_cent_model = block['SigmaCentY1CUDAIntegrand', "vals"]
    NC_miscent_model = block['NCMiscentY1CUDAIntegrand', "vals"]
    NC_cent_model = block['NCCentY1CUDAIntegrand', "vals"]
    miscent_fraction = block['cluster_abundance', 'fcen']

    # TODO: Consider making the value of too_small a configuration parameter
    # of the module.
    too_small = 0.1
    if (any_model_too_small(too_small, profiles_miscent_model, profiles_cent_model, NC_miscent_model, NC_cent_model)):
        block[section_names.likelihoods, "SigmaMort_Like_like"] = -np.inf
        return 0

    # Recover what we created in setup.
    data_vector, covmat, ncs, nc_covmat, rad, Redges = config
   
    sigma = profiles_miscent_model/NC_miscent_model * miscent_fraction + (profiles_cent_model/NC_cent_model)*(1-miscent_fraction)
    model_ncs = NC_miscent_model * miscent_fraction + (NC_cent_model)*(1-miscent_fraction)
    model_ncs = model_ncs[:, 0]
    ds, aver_ds = convert_s2ds(rad, sigma, Redges)
    model_vector = aver_ds.flatten()

    delta = data_vector - model_vector
    weight = np.linalg.inv(covmat)
    loglike1 = -0.5 * np.dot(delta, np.dot(weight, delta))
    
    delta2 = ncs - model_ncs
    weight2 = np.linalg.inv(nc_covmat)
    loglike2 = -0.5 * np.dot(delta2, np.dot(weight2, delta2))
    loglike = loglike1 + loglike2

    block[section_names.likelihoods, 'SigmaMort_Like_like'] = loglike

    #print(block['cluster_mass_profile', "concentration"], block['cosmological_parameters', "omega_m"], block['cosmological_parameters', "sigma_8"], loglike1, loglike2, loglike, delta, delta2, model_ncs, ncs, model_vector)
    #plt.plot(model_vector, 'k--');plt.plot(data_vector, 'b:');plt.show()
    #plt.plot(ncs, 'k--');plt.plot(model_ncs, 'b:');plt.show()
    #make_plots(3, 4, 18, model_vector, data_vector, covmat)
    #make_plots_compare(3, 4, 18, model_vector, data_vector, covmat)
    return 0

def cleanup(config):
    #nothing to clean up
    return 0

def any_model_too_small(small, *models):
    """Take any number of np.ndarray objects, and check that all entries are
    'not too small'.

    Return True if any entry is too small, and False otherwise.
    """
    for model in models:
        if np.any(model < small):
            return True
    return False

def sig_gamma(rad, profiles, xx0):
    ds= np.zeros([profiles.shape[0], len(rad)])
    #print(ds.shape)
    for jj in range(profiles.shape[0]):
        for ii in range(len(rad)):
            ds[jj, ii] = np.trapz(profiles[jj, :ii]*rad[jj, :ii], rad[jj, :ii])*2.0/(rad[ii])**2 - profiles[jj, ii]
    f=interp1d(rad, ds)
    mean=ds -f(xx0)*xx0**2/rad**2
    return mean

def convert_s2ds(rad, profiles, Redges):
    ds= np.zeros([profiles.shape[0], len(rad)])
    #print(ds.shape, profiles.shape)
    for jj in range(profiles.shape[0]):
        for ii in range(len(rad)):
            ds[jj, ii] = np.trapz(profiles[jj, :ii]*rad[:ii], rad[:ii])*2.0/(rad[ii])**2 - profiles[jj, ii]


    res = np.zeros([profiles.shape[0], len(Redges)-1])
    # MFP: Commented out this double loop, because it appears to have no purpose except to print -- and we don't
    # want the printout in production running.
    #for jj in range(profiles.shape[0]):
    #    prof_interp=interp1d(rad, ds[jj, :])
    #    prof_func=lambda x: prof_interp(x) * x
    #    for ii in range(len(Redges)-1):
    #        print(jj, ii, integrate.quad(prof_func, Redges[ii], Redges[ii+1])*2/((Redges[ii+1]**2) - (Redges[ii]**2)))
    return ds, res


def assemble_vector(profiles_model, NCs_model, radii, snapshot_zs):
    averaged_profiles = profiles_model/NCs_model
    for ii in range(len(snapshot_zs)):
        if ii == 0:
           model_vec = NCs_model[:, 0]
        else:
           model_vec=np.append(model_vec, NCs_model[:,len(radii)*ii])
    model_vec=np.append(model_vec, averaged_profiles)

    return model_vec

def make_plots(ngriddim, nvoldim, nrad, model_vector, data_vector, covmat):
    ngriddim = int(ngriddim)
    nvoldim = int(nvoldim)
    nrad = int(nrad)
    unc = np.diag( np.sqrt(covmat) )
    fig, axs = plt.subplots(ngriddim, nvoldim+1, figsize=[12, 8])

    for j in range(ngriddim):
        for i in range(nvoldim):
            ndim1=j*nvoldim
            ndim2=(j+1)*nvoldim
            axs[j, nvoldim].plot(model_vector[ndim1:ndim2])
            axs[j, nvoldim].errorbar(range(nvoldim), data_vector[ndim1:ndim2], yerr=unc[ndim1:ndim2], fmt='ko')

            dim1 = i*ngriddim*nrad + j*nrad + ngriddim*nvoldim
            dim2 = i*ngriddim*nrad + (j+1)*nrad + ngriddim*nvoldim
            axs[j, i].plot(range(nrad), model_vector[dim1:dim2], 'r--')
            axs[j, i].errorbar(range(nrad), data_vector[dim1:dim2], yerr=unc[dim1:dim2], fmt='ko')
            axs[j, i].set_yscale('log')
    plt.show()


def make_plots_compare(ngriddim, nvoldim, nrad, model_vector, data_vector, covmat):
    ngriddim = int(ngriddim)
    nvoldim = int(nvoldim)
    nrad = int(nrad)
    unc = np.diag( np.sqrt(covmat) )
    fig, axs = plt.subplots(ngriddim, nvoldim+1, figsize=[12, 8])

    for j in range(ngriddim):
        for i in range(nvoldim):
            ndim1=j*nvoldim
            ndim2=(j+1)*nvoldim
            axs[j, nvoldim].errorbar(range(nvoldim), data_vector[ndim1:ndim2]/model_vector[ndim1:ndim2], yerr=unc[ndim1:ndim2]/model_vector[ndim1:ndim2], fmt='ko')

            dim1 = i*ngriddim*nrad + j*nrad + ngriddim*nvoldim
            dim2 = i*ngriddim*nrad + (j+1)*nrad + ngriddim*nvoldim
            axs[j, i].errorbar(range(nrad), data_vector[dim1:dim2]/model_vector[dim1:dim2], yerr=unc[dim1:dim2]/model_vector[dim1:dim2], fmt='ko')
            axs[j, i].set_ylim([0.5, 2.0])
            axs[j, nvoldim].set_ylim([0.5, 2.0])
    plt.show()





