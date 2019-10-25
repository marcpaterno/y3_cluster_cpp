import numpy as np
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section
import matplotlib.pyplot as plt

def setup(options):
    section = option_section

    #read in the profile_file, need to make this real
    profile_file = options[section,"profile_file"]
    profiles=np.genfromtxt(profile_file)

    # read in the NC file, need to make this real
    NC_file = options[section,"NC_file"]
    NCs=np.genfromtxt(NC_file)
    data_vector=assemble_vector(profiles*NCs, NCs)
    print(NCs.shape, data_vector.shape, profiles.shape)
    #read in the profile_errfile, need to make this real

    profile_file_err = options[section,"profile_err_file"]
    profiles_err=np.genfromtxt(profile_file_err)

    # read in the NCerr file
    NC_file_err = options[section,"NC_err_file"]
    NCs_err=np.genfromtxt(NC_file_err)
    cov_vector=assemble_vector(profiles_err*NCs_err, NCs_err)
    #covariance_file= options[section,"covariance_file"]
    #Need to make this real
    covmat = np.diag(cov_vector**2)
     
    print(NCs, NCs_err, cov_vector, profiles_err, covmat)

    return data_vector, covmat

def execute(block, config):
    data_vector, covmat = config

    #read in the model values at this sample point. 
    profile_blockname ='snapshotscalarintegrand' 
    profiles_model = block[profile_blockname, "vals"]
    NCs_blockname ='snapshotscalarncintegrand' 
    NCs_model = block[NCs_blockname, "vals"]
    model_vector = assemble_vector(profiles_model, NCs_model)

    delta = data_vector - model_vector
    weight = np.linalg.inv(covmat)
    loglike = -0.5 * np.dot(delta, np.dot(weight, delta))
    block[section_names.likelihoods, 'SnapshotScalarLike_like'] = loglike

    print(block['cluster_mass_profile', "concentration"], block['cosmological_parameters', "omega_m"], block['cosmological_parameters', "sigma_8"], loglike)
    '''plt.plot( model_vector );
    plt.yscale('log')
    plt.errorbar(range(len(model_vector)), data_vector, yerr=np.sqrt(np.diag(covmat)), fmt='r.')
    plt.show()'''
    return 0

def cleanup(config):
    #nothing to clean up
    return 0

def assemble_vector(profiles_model, NCs_model):
    averaged_profiles = profiles_model/NCs_model
    model_vec=np.array([NCs_model[:,0], NCs_model[:,19], NCs_model[:,38]])
    model_vec=np.append(model_vec, averaged_profiles)

    return model_vec
