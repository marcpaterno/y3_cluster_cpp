import numpy as np
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
import matplotlib.pyplot as plt

def setup(options):
    section = option_section

    #read in prediction names
    profile_blockname = options[section,"profile_blockname"]
    NCs_blockname = options[section,"NCs_blockname"]
    snapshot_zs = options[section, "snapshot_zs"]
    snapshot_lnm = options[section, "lnm_low"]
    snapshot_radii = options[section, "radii"]
    #read in the profile_file, need to make this real
    profile_file = options[section,"profile_file"]
    profiles=np.genfromtxt(profile_file)

    # read in the NC file, need to make this real
    NC_file = options[section,"NC_file"]
    NCs=np.genfromtxt(NC_file)
    data_vector=assemble_vector(profiles*NCs, NCs, snapshot_radii, snapshot_zs)
    print(NCs.shape, data_vector.shape, profiles.shape)
    #read in the profile_errfile, need to make this real

    profile_file_err = options[section,"profile_err_file"]
    profiles_err=np.genfromtxt(profile_file_err)

    # read in the NCerr file
    NC_file_err = options[section,"NC_err_file"]
    NCs_err=np.genfromtxt(NC_file_err)
    cov_vector=assemble_vector(profiles_err*NCs_err, NCs_err, snapshot_radii, snapshot_zs)
    #covariance_file= options[section,"covariance_file"]
    #Need to make this real
    covmat = np.diag(cov_vector**2)


    return data_vector, covmat, profile_blockname, NCs_blockname, snapshot_zs, snapshot_lnm, snapshot_radii

def execute(block, config):
    data_vector, covmat, profile_blockname, NCs_blockname, snapshot_zs, snapshot_lnm, snapshot_radii = config

    #read in the model values at this sample point. 
    profiles_model = block[profile_blockname, "vals"]
    NCs_model = block[NCs_blockname, "vals"]
    model_vector = assemble_vector(profiles_model, NCs_model, snapshot_radii, snapshot_zs)

    delta = data_vector - model_vector
    weight = np.linalg.inv(covmat)
    loglike = -0.5 * np.dot(delta, np.dot(weight, delta))
    block[section_names.likelihoods, 'SnapshotScalarLike_like'] = loglike

    print(block['cluster_mass_profile', "concentration"], block['cosmological_parameters', "omega_m"], block['cosmological_parameters', "sigma_8"], loglike)
    make_plots(3, 4, 18, model_vector, data_vector, covmat)
    make_plots_compare(3, 4, 18, model_vector, data_vector, covmat)
    return 0

def cleanup(config):
    #nothing to clean up
    return 0

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





