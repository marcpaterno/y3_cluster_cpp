# Apply Boost Factor Correction to Shear Signal
# Author: Arwa's code with Claude's modifications
# Reference: McClintock et al. 2019 (arXiv:1805.00039)

# This module divides the theoretical shear by B(R) before the likelihood
# comparison, so the model prediction matches what we'd observe with
# source contamination.

import numpy as np

# CosmoSIS imports - only needed when running within CosmoSIS
try:
    from cosmosis.datablock import option_section, names
    from utils import find_cosmosis_sections
    HAS_COSMOSIS = True
except ImportError:
    HAS_COSMOSIS = False
    option_section = None
    names = None

# Import bin configuration (only needed within CosmoSIS)
try:
    from setup_bins import zmeans_ij, lambda_means_ij
except ImportError:
    zmeans_ij = None
    lambda_means_ij = None

if HAS_COSMOSIS:
    cosmo = names.cosmological_parameters
else:
    cosmo = None


def boost_factor_model(R, rs, b0):
    """
    Compute boost factor B(R) using NFW-like profile.

    Parameters
    ----------
    R : array_like
        Radial distances (Mpc/h)
    rs : float
        Scale radius parameter (Mpc/h)
    b0 : float
        Amplitude parameter

    Returns
    -------
    B : ndarray
        Boost factor values B(R) >= 1
    """
    x = np.atleast_1d(R / rs).astype(float)
    B = np.zeros_like(x, dtype=float)

    tol = 1e-6

    mask_near1 = np.abs(x - 1) < tol
    mask_gt1 = (x > 1) & ~mask_near1
    mask_lt1 = (x < 1) & ~mask_near1

    # x > 1: arctan regime
    if np.any(mask_gt1):
        x_gt1 = x[mask_gt1]
        sqrt_term = np.sqrt(x_gt1**2 - 1)
        fx = np.arctan(sqrt_term) / sqrt_term
        B[mask_gt1] = 1 + b0 * (1 - fx) / (x_gt1**2 - 1)

    # x < 1: arctanh regime
    if np.any(mask_lt1):
        x_lt1 = x[mask_lt1]
        sqrt_term = np.sqrt(1 - x_lt1**2)
        fx = np.arctanh(sqrt_term) / sqrt_term
        B[mask_lt1] = 1 + b0 * (1 - fx) / (x_lt1**2 - 1)

    # x ~ 1: analytic limit
    if np.any(mask_near1):
        B[mask_near1] = (b0 + 3) / 3

    # Handle NaN/inf
    B = np.where(np.isnan(B) | np.isinf(B), (b0 + 3) / 3, B)

    return B


def setup(options):
    """
    Setup function - read configuration from .ini file.

    Options:
    - n_lambda_bins: Number of richness bins (default: 4)
    - n_z_bins: Number of redshift bins (default: 3)
    - shear_section: Datablock section for shear (default: auto-detect)
    """
    section = option_section

    n_lambda_bins = options.get_int(section, "n_lambda_bins", default=4)
    n_z_bins = options.get_int(section, "n_z_bins", default=3)

    # Output section name for corrected shear
    output_section = options.get_string(section, "output_section", default="shear_boosted")

    config = {
        'n_lambda_bins': n_lambda_bins,
        'n_z_bins': n_z_bins,
        'output_section': output_section,
    }

    print(f"ApplyBoostFactor: {n_lambda_bins} lambda bins x {n_z_bins} z bins")
    return config


def execute(block, config):
    """
    Apply boost factor correction to shear signal.

    Reads shear from datablock, computes B(R) for each bin,
    divides shear by B(R), writes corrected shear.
    """
    n_lambda_bins = config['n_lambda_bins']
    n_z_bins = config['n_z_bins']
    output_section = config['output_section']

    # Find the shear section in datablock
    shear_section = find_cosmosis_sections(block, "shear")[0]

    # Read shear data
    shear_cen = block[shear_section, "shear_cen"]  # shape: (Nlbins * Nzbins, Nrbins)
    r_shear = block[shear_section, "r"]  # radii in Mpc/h

    Nlbins_Nzbins, Nrbins = shear_cen.shape

    # Verify dimensions
    expected_bins = n_lambda_bins * n_z_bins
    if Nlbins_Nzbins != expected_bins:
        print(f"Warning: shear has {Nlbins_Nzbins} bins, expected {expected_bins}")

    # Apply boost correction to each bin
    shear_corrected = np.zeros_like(shear_cen)

    for ij in range(Nlbins_Nzbins):
        # Determine lambda and z bin indices
        # Assuming ordering: (l0,z0), (l0,z1), (l0,z2), (l1,z0), ...
        l_bin = ij // n_z_bins
        z_bin = ij % n_z_bins

        # Read boost parameters for this bin
        param_suffix = f"l{l_bin}_z{z_bin}"

        try:
            logrs = block["boost_factor_params", f"logrs_{param_suffix}"]
            logb0 = block["boost_factor_params", f"logb0_{param_suffix}"]
        except:
            # If parameters not found, use default (no boost correction)
            print(f"Warning: boost params not found for {param_suffix}, using B=1")
            shear_corrected[ij] = shear_cen[ij]
            continue

        rs = 10**logrs
        b0 = 10**logb0

        # Compute boost factor B(R)
        B_R = boost_factor_model(r_shear, rs, b0)

        # Apply correction: gamma_model_effective = gamma_model / B(R)
        shear_corrected[ij] = shear_cen[ij] / B_R

    # Write corrected shear to datablock
    # Overwrite the original shear section so likelihood uses corrected values
    block[shear_section, "shear_cen"] = shear_corrected

    # Also store in separate section for diagnostics
    block[output_section, "shear_cen"] = shear_corrected
    block[output_section, "shear_uncorrected"] = shear_cen
    block[output_section, "r"] = r_shear

    return 0


def cleanup(config):
    pass
