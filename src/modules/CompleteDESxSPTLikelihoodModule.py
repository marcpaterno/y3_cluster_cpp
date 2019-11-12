import numpy as np
from cosmosis.datablock import names, option_section


def setup(options):
    """Nothing to read in or setup to perform"""
    return 0


def execute(block, config):
    r"""Compute the log-likelihood value and save it to the datablock

    \ln\mathcal{L} = \sum_{i=1}^N \ln P(\lambda^{obs}, \xi, z^{obs})
        - \int_{\lambda_N}^\infty d\lambda \int_{survey} dz^{obs} P(\lambda^{obs}, z^{obs})
    """
    # Load the precompute integration results
    # abundances = block.get_double_array_1d("abundance_integral", "vals")
    abundances = block["abundance_integral", "vals"]
    selection_term = block["selection_integral", "vals"]

    # Compute the total lnlikelihood and finish up
    lnlike = np.sum(np.log(abundances)) - selection_term
    block.put_double("richest_clusters_likelihood", "lnlike", lnlike)
    return 0


def cleanup(config):
    pass
