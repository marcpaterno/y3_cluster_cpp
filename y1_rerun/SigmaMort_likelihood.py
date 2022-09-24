import numpy as np
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.interpolate import interp1d

from mpi4py import MPI

def setup(options):
    section = option_section
    # data vector 
    profiles = np.genfromtxt('data_files/wl_data_vector.txt')
    profiles_err = np.genfromtxt('data_files/wl_data_vector_err.txt')
    covmat = np.genfromtxt('data_files/wl_cov.txt')
    ncs = np.genfromtxt('data_files/nc.txt', delimiter=',')
    nc_covmat = np.genfromtxt('data_files/Cov_ij_bestfit_DESY1_105.txt')
    #ncs = ncs[0, :]
    nc_covmat = nc_covmat[:12, :12]

    # model module     
    '''profile_blockname = options[section,"profile_blockname"]
    NCs_blockname = options[section,"NCs_blockname"]
    snapshot_zs = options[section, "snapshot_zs"]
    snapshot_lnm = options[section, "lnm_low"]
    snapshot_radii = options[section, "radii"]'''
    rad = options[section, "radius"]
    Redges = np.logspace(np.log10(0.0323), np.log10(30.), num=15+1)
    return profiles, covmat, ncs, nc_covmat, rad, Redges, profiles_err, 

def execute(block, config):
    print(f"In SigmaMort_likelihood.execute: mpi rank is {MPI.COMM_WORLD.rank} total number of ranks is: {MPI.COMM_WORLD.size}")
    data_array, covmat, ncs, nc_covmat, rad, Redges, data_vector_err = config

    #read in the model values at this sample point.
    # We look for CUDA algorithm results first; if that fails, we look for CPU-based
    # algorithm results. If that fails the module fails. This implementation assumes
    # that we are either using all CUDA algorithms or all CPU algorithms. If that is
    # not sufficient, we need a slightly more complex bit of logic here.
    try:
        profiles_miscent_model = block['SigmaMiscentY1MortCUDAScalarIntegrand', "vals"]
        profiles_cent_model = block['SigmaCentY1MortCUDAScalarIntegrand', "vals"]
        NC_miscent_model = block['NCMiscentY1MortCUDAScalarIntegrand', "vals"]
        NC_cent_model = block['NCCentY1MortCUDAScalarIntegrand', "vals"]
    except:
        profiles_miscent_model = block['SigmaMiscentY1MortScalarIntegrand', "vals"]
        profiles_cent_model = block['SigmaCentY1MortScalarIntegrand', "vals"]
        NC_miscent_model = block['NCMiscentY1MortScalarIntegrand', "vals"]
        NC_cent_model = block['NCCentY1MortScalarIntegrand', "vals"]

    miscent_fraction = block['cluster_abundance', 'fcen'] 
    
    '''
    profiles_miscent_model = np.genfromtxt('sigma_y1_output/sigmamiscenty1mortscalarintegrand/vals.txt')
    profiles_cent_model = np.genfromtxt('sigma_y1_output/sigmacenty1mortScalarintegrand/vals.txt')
    NC_miscent_model = np.genfromtxt('sigma_y1_output/ncmiscentY1mortscalarintegrand/vals.txt')
    NC_cent_model = np.genfromtxt('sigma_y1_output/nccenty1mortscalarintegrand/vals.txt')
    miscent_fraction = 0.7 #np.genfromtxt('sigma_y1_output/cluster_abundance/fcen.txt') 
    '''

    sigma = profiles_miscent_model/NC_miscent_model * miscent_fraction + (profiles_cent_model/NC_cent_model)*(1-miscent_fraction)
    ncs_vector = NC_miscent_model * miscent_fraction + (NC_cent_model)*(1-miscent_fraction)
    #sigma = profiles_miscent_model/NC_miscent_model # * miscent_fraction # + (profiles_cent_model/NC_cent_model)*(1-miscent_fraction)
    #ncs_vector = NC_miscent_model # * miscent_fraction # + (NC_cent_model)*(1-miscent_fraction)
    #sigma = (profiles_cent_model/NC_cent_model)
    #ncs_vector = (NC_cent_model)
    model_ncs = np.array([])
    ngrid = sigma.shape[1]/len(rad)
    ngrid = int(ngrid)
    for kk in range(ngrid):
        model_ncs = np.append(model_ncs, ncs_vector[:, kk*(len(rad))])
    ds, aver_ds = convert_s2ds(rad, sigma, Redges)
    model_vector = aver_ds.flatten()
    
    data_vector = data_array.flatten()
    delta = data_vector - model_vector
    weight = np.linalg.inv(covmat)
    loglike1 = -0.5 * np.dot(delta, np.dot(weight, delta))
    
    loglike1 = 0.0
    print(ncs.shape, model_ncs.shape, nc_covmat.shape)
    ncs = ncs.flatten()
    delta2 = ncs - model_ncs
    weight2 = np.linalg.inv(nc_covmat)
    loglike2 = -0.5 * np.dot(delta2, np.dot(weight2, delta2))
    loglike = loglike1 + loglike2
    block[section_names.likelihoods, 'SigmaMort_Like_like'] = loglike

    #print(block['cluster_mass_profile', "concentration"], block['cosmological_parameters', "omega_m"], block['cosmological_parameters', "sigma_8"], loglike1, loglike2, loglike, delta, delta2, model_ncs, ncs, model_vector)
    #plt.plot(model_vector, 'k--');plt.plot(data_vector, 'b:');plt.show()
    #plt.plot(ncs, 'k--');plt.plot(model_ncs, 'b:');plt.show()
    #make_plots(3, 4, Redges, aver_ds, data_array,  data_vector_err, ncs, model_ncs, nc_covmat)
    #make_plots_compare(3, 4, 18, model_vector, data_vector, covmat)

    return 0

def cleanup(config):
    #nothing to clean up
    return 0

def sig_gamma(rad, profiles, xx0):
    ds= np.zeros([profiles.shape[0], len(rad)])
    print(ds.shape)
    for jj in range(profiles.shape[0]):
        for ii in range(len(rad)):
            ds[jj, ii] = np.trapz(profiles[jj, :ii]*rad[jj, :ii], rad[jj, :ii])*2.0/(rad[ii])**2 - profiles[jj, ii]
    f=interp1d(rad, ds)
    mean=ds -f(xx0)*xx0**2/rad**2
    return mean

def convert_s2ds(rad, profiles, Redges):
    #print(ds.shape, profiles.shape)
    ngrid=int(profiles.shape[1]/len(rad))
    ds= np.zeros([profiles.shape[0]*ngrid, len(rad)])
    for jj in range(profiles.shape[0]):
        for kk in range(ngrid):
            for ii in range(len(rad)):
                profiles_ind = kk * len(rad)
                ds_int = np.trapz(profiles[jj, profiles_ind:(profiles_ind+ii+1)]*rad[:(ii+1)], rad[:(ii+1)])*2.0/(rad[ii])**2 - profiles[jj, (profiles_ind+ii)]
                ds[kk*profiles.shape[0]+jj, ii] = ds_int


    res = np.zeros([ds.shape[0], len(Redges)-1])
    for jj in range(ds.shape[0]):
        prof_interp=interp1d(rad, ds[jj, :])
        #prof_interp=interp1d(rad, profiles[jj, :])
        prof_func=lambda x: prof_interp(x) * x
        for ii in range(len(Redges)-1):
            # averaging DS ranges
            res_int, unc_int = integrate.quad(prof_func, Redges[ii], Redges[ii+1])
            res_int = res_int *2/((Redges[ii+1]**2) - (Redges[ii]**2))

            #print(jj, ii, res_int)
            res[jj, ii] = res_int
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

def make_plots(ngriddim, nvoldim, rad_edges, model_vector, data_vector, data_err, ncs, model_ncs, nc_covmat):

    rad = 0.5*(rad_edges[1:]**2 - rad_edges[:-1]**2)/(rad_edges[1:] - rad_edges[:-1])
    ngriddim = int(ngriddim)
    nvoldim = int(nvoldim)
    nrad = int(len(rad))
    unc = data_err
    unc_nc = np.diag( np.sqrt(nc_covmat) )
    fig, axs = plt.subplots(2, ngriddim, figsize=[12, 8])

    zbin_labels = [r"Redshift Bin [0.2;0.35)", r"Redshift Bin [0.35;0.5)", r"Redshift Bin [0.5;0.65)"]
    lam_labels = [r"$\lambda \in [20;30)$", r"$\lambda \in [30;45)$", r"$\lambda \in [45;60)$", r"$\lambda \in [60;\infty)$"]
    for j in range(ngriddim):
        axs[0, j].errorbar(range(nvoldim), model_ncs[j*nvoldim:j*nvoldim+nvoldim], fmt='s', color = 'tab:pink', label = 'Models')
        axs[0, j].errorbar(range(nvoldim), ncs[j*nvoldim:j*nvoldim+nvoldim], yerr=unc_nc[j*nvoldim:j*nvoldim+nvoldim], fmt='o', color = 'gray', label = 'Data')
        axs[0, j].set_xticks(range(nvoldim))
        axs[0, j].set_xticklabels(lam_labels) 
        for i in range(nvoldim):
            #ndim1=j*nvoldim
            #ndim2=(j+1)*nvoldim
            #axs[j, nvoldim].plot(model_vector[ndim1:ndim2])
            #print(model_ncs, ncs, unc_nc)
            axs[1, j].plot(rad, model_vector[j*nvoldim+i], label = "Models "+lam_labels[i])
            #print(dim1, dim2, len(unc), model_vector[dim1:dim2])
            dim1 = i * len(rad)
            dim2 = (i+1) * len(rad)
            axs[1, j].errorbar(rad, data_vector[j*nvoldim+i], yerr=unc[j*nvoldim+i], fmt='o', color = 'gray')
            axs[1, j].set_yscale('log')
            axs[1, j].set_xscale('log')
            axs[1, j].set_xlim([0.2, 20])
            axs[1, j].set_ylim([2, 2000])
            
            axs[0, j].set_yscale('log')
            axs[1, j].set_xlabel(r'$R$ [Comoving] Mpc/$h$')
        axs[0, j].set_title(zbin_labels[j])
    axs[0, 0].set_ylabel('Number Counts')
    axs[1, 0].set_ylabel(r'$\Delta \Sigma$')
    axs[1, 2].legend(loc = 0)
    axs[0, 2].legend(loc = 0)
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





