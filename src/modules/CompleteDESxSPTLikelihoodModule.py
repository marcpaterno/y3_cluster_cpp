import numpy as np
from cosmosis.datablock import names, option_section


def setup(options):
    """Nothing to read in or setup to perform"""
    return 0


def execute(block, config):
    r"""Compute the log-likelihood value and save it to the datablock

    \ln\mathcal{L} = \sum_{i=1}^N \ln P(\lambda^{obs}, \xi, z^{obs})
        - \int_{\lambda_N}^\infty d\lambda \int_{survey} dz^{obs} P(\lambda^{obs}, z^{obs})

    Approximations
    1. \xi = \zeta (no maximization bias)
    2. ztrue integral is done over P(zo|zt)A(zt)dv/dztdo only
    3. P(lo|lt,z), n(M,z), P(lt,zeta|M,z) evaluated at zobs
    """
    #########################################
    # Compute the first (abundance) term
    #########################################
    # Load the ztrue integration
    ztrue_status = block["ztrue_integration", "status"]
    ztrue_integral = block["ztrue_integration", "vals"]
    if np.any(ztrue_status > 0.9):
        raise RuntimeError("ztrue integration unconverged.")

    # Load the lamtrue and lnM integration
    lamtrue_lnm_integral = block["abundance_integral", "vals"]

    # Abundance integral - Need to revise after P(xi|zeta) added in
    abundances = ztrue_integral * lamtrue_lnm_integral


    #########################################
    # Load the second (selection) term
    #########################################
    selection_term = block["selection_integral", "vals"]

    #########################################
    # Compute the total lnlikelihood
    #########################################
    lnlike = np.sum(np.log(abundances)) - selection_term
    block.put_double("richest_clusters_likelihood", "lnlike", lnlike)
    return 0


def cleanup(config):
    pass
