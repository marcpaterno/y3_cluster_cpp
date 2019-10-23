import numpy as np
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section


def setup(options):
    section = option_section

    #read in the profile_file, need to make this real
    profile_file = options[section,"profile_file"]
    profiles=np.genfromtxt(profile_file)

    # read in the NC file, need to make this real
    NC_file = options[section,"NC_file"]
    NCs=np.genfromtxt(NC_file)
    data_vector=assemble_vector(profiles*NCs, NCs)

    #covariance_file= options[section,"covariance_file"]
    #Need to make this real
    covmat = np.eye(len(data_vector)) 
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
    loglike = np.dot(delta, np.dot(weight, delta))
    block[section_names.likelihoods, 'SnapshotScalarLike_like'] = loglike

    return 0

def cleanup(config):
    #nothing to clean up
    return 0

def assemble_vector(profiles_model, NCs_model):
    averaged_profiles = profiles_model/NCs_model
    model_vec=np.append(NCs_model[0, 0], NCs_model[1, 0])
    model_vec=np.append(model_vec, profiles_model)

    return model_vec
