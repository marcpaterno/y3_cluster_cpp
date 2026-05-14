"""Gaussian likelihood for the mock_mcmc_cp_camb.ini production pipeline.

Compares fiducial (data, invcov) against two observables from the
Costanzi-2026 pipeline:

    Number counts         numcountssel/vals          -> 12
    Tangential shear      gamma_t^theory(R)          -> 120
        = <gamma_t^1h>(R) + gamma_t^prj(R)
        = (shear1hsel/vals / numcountssel/vals) + shear_prj/vals

The two shear contributions are added linearly — which is valid only
because we dropped the reduced-shear 1/(1 - sigma sci) denominator
(see src/modules/num_counts_sel/lensing_weights.hh and
src/models/sigma_prj_t.hh, both updated 2026-05-11).

shear1hsel/vals is the N_i-weighted integral of gamma_t^1h
(N_i[g]/N_i is the per-bin average); divide by numcountssel/vals to
get the bin-averaged 1-halo shear.  shear_prj/vals is already a
per-(lambda_bin, zob) shear ready to add in.

Diagonal invcov is the default (1-D array per observable); a dense
invcov (2-D array) is also accepted and used as a matmul.

logL = -0.5 * sum_obs delta^T C^-1 delta, summed over NC and Shear.

Writes block["likelihoods", "likelihoods_like"].
"""

from __future__ import annotations

import numpy as np
from cosmosis.datablock import option_section

_NC_N_BINS = 12
_SHEAR_N_R = 10         # radii per bin
_SHEAR_N = _NC_N_BINS * _SHEAR_N_R   # 120

OBS = [
    ("NC",    _NC_N_BINS),
    ("Shear", _SHEAR_N),
]


def _chi2(delta: np.ndarray, invcov: np.ndarray) -> float:
    if invcov.ndim == 1:
        return float(np.sum(delta * delta * invcov))
    if invcov.ndim == 2:
        return float(delta @ invcov @ delta)
    raise ValueError(f"invcov ndim must be 1 or 2, got {invcov.ndim}")


def setup(options):
    fname = options[option_section, "filename"]
    vec = np.load(fname, allow_pickle=False)
    config = {"filename": fname,
              "verbose": bool(options.get_bool(option_section, "verbose",
                                                default=False))}
    for name, expected_n in OBS:
        d = np.asarray(vec[f"data_{name}"]).ravel()
        ic = np.asarray(vec[f"invcov_{name}"])
        if d.size != expected_n:
            raise ValueError(
                f"likelihood_cp: data_{name} has size {d.size}, "
                f"expected {expected_n}")
        if ic.ndim == 1 and ic.size != expected_n:
            raise ValueError(
                f"likelihood_cp: invcov_{name} diagonal has size {ic.size}, "
                f"expected {expected_n}")
        if ic.ndim == 2 and ic.shape != (expected_n, expected_n):
            raise ValueError(
                f"likelihood_cp: invcov_{name} dense has shape {ic.shape}, "
                f"expected ({expected_n},{expected_n})")
        config[f"data_{name}"] = d
        config[f"invcov_{name}"] = ic
    print(f"[likelihood_cp] loaded mock DV from {fname}: "
          f"NC={config['data_NC'].size}, Shear={config['data_Shear'].size}")
    return config


def _shear_theory(block) -> np.ndarray:
    """Build the summed tangential-shear theory vector (length 120).

    gamma_t^theory(R | i,j) = <gamma_t^1h>_i(R) + gamma_t^prj(R | i,j)

    shear1hsel/vals is N_i-weighted (shape 120); divide entry-wise by
    the NumCountsSel integral (shape 12, broadcast across the 10 R
    points per bin) to get the per-cluster average, then add the
    projection piece.
    """
    NC = np.asarray(block["numcountssel", "vals"]).ravel()
    S1h_Ni = np.asarray(block["shear1hmissel", "vals"]).ravel()
    Sprj = np.asarray(block["shear_prj", "vals"]).ravel()
    if NC.size != _NC_N_BINS:
        raise ValueError(
            f"likelihood_cp: numcountssel/vals size {NC.size} != "
            f"{_NC_N_BINS}")
    if S1h_Ni.size != _SHEAR_N:
        raise ValueError(
            f"likelihood_cp: shear1hmissel/vals size {S1h_Ni.size} != "
            f"{_SHEAR_N}")
    if Sprj.size != _SHEAR_N:
        raise ValueError(
            f"likelihood_cp: shear_prj/vals size {Sprj.size} != "
            f"{_SHEAR_N}")
    # shear1hsel wall = (bin_index fast, r_perp slow) for Cartesian
    # product; NumCountsSel wall is just bin_index.  The two module
    # outputs share the same 12-bin ordering, so we can tile NC across
    # the 10 R-per-bin axis.
    NC_tile = np.repeat(NC, _SHEAR_N_R)
    bad = NC_tile <= 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        S1h_avg = np.where(bad, 0.0, S1h_Ni / NC_tile)
    return S1h_avg + Sprj


def execute(block, config):
    logL = 0.0
    parts = {}

    # NumCounts — direct Gaussian on the 12-bin vector.
    NC_theory = np.asarray(block["numcountssel", "vals"]).ravel()
    delta_NC = config["data_NC"] - NC_theory
    parts["NC"] = -0.5 * _chi2(delta_NC, config["invcov_NC"])
    logL += parts["NC"]

    # Shear — theory = <gamma_t^1h> + gamma_t^prj (length 120).
    Shear_theory = _shear_theory(block)
    delta_Shear = config["data_Shear"] - Shear_theory
    parts["Shear"] = -0.5 * _chi2(delta_Shear, config["invcov_Shear"])
    logL += parts["Shear"]

    if config["verbose"]:
        print(f"[likelihood_cp] logL={logL:.4e}  "
              f"(NC={parts['NC']:.3e}, Shear={parts['Shear']:.3e})")

    block["likelihoods", "likelihoods_like"] = float(logL)
    return 0


def cleanup(config):
    return 0
