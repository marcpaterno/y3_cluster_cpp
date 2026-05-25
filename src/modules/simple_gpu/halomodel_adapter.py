"""Adapter module: copies haloModel section to deltasigma section.

The GPU modules (SIG_MAX, GAMMA_MAX) expect data in section 'deltasigma'
with specific key names. This adapter copies from 'haloModel' section
(written by buildWpGammat.py) to 'deltasigma' with the expected names.

Mapping:
  haloModel/Sigma_nfw  -> deltasigma/sigma_1
  haloModel/Sigma_hh   -> deltasigma/sigma_2
  haloModel/dSigma_nfw -> deltasigma/dsigma_1
  haloModel/dSigma_hh  -> deltasigma/dsigma_2
  haloModel/bias       -> deltasigma/bias
  haloModel/r_sigma    -> deltasigma/r_sigma_deltasigma
  haloModel/lnM        -> deltasigma/lnM
  haloModel/z          -> deltasigma/z
"""

import numpy as np


def setup(options):
    """No configuration needed."""
    return {}


def execute(block, config):
    """Copy haloModel section to deltasigma section with renamed keys."""

    src = "haloModel"
    dst = "deltasigma"

    # Copy with key name mapping
    # Note: sigma_1 needs to be transposed for SIG_MAX (expects r x lnM)
    # buildWpGammat outputs Sigma_nfw as (M_bins x R_bins)

    sigma_nfw = block[src, "Sigma_nfw"]
    sigma_hh = block[src, "Sigma_hh"]
    dsigma_nfw = block[src, "dSigma_nfw"]
    dsigma_hh = block[src, "dSigma_hh"]
    bias = block[src, "bias"]
    r_sigma = block[src, "r_sigma"]
    lnM = block[src, "lnM"]
    z = block[src, "z"]

    # Transpose sigma arrays: (M, R) -> (R, M) for SIG_MAX interpolator
    # SIG_MAX expects sigma_1(r, lnM) with r as first axis
    block[dst, "sigma_1"] = sigma_nfw.T
    block[dst, "sigma_2"] = sigma_hh.T
    block[dst, "dsigma_1"] = dsigma_nfw.T
    block[dst, "dsigma_2"] = dsigma_hh.T

    # Bias is (M, z), SIG_MAX expects (z, lnM)
    block[dst, "bias"] = bias.T

    # Copy axes
    block[dst, "r_sigma_deltasigma"] = r_sigma
    block[dst, "lnM"] = lnM
    block[dst, "z"] = z

    return 0


def cleanup(config):
    return 0
