"""Combine outputs of the four NOperator[F] modules into effective per-bin
observables:

    M_eff[i]          = N_i[M]           / N_i[1]
    b_eff[i]          = N_i[b(M,z)]      / N_i[1]
    DeltaSigma_1h[i, r] = N_i[ΔΣ_nfw(r)] / N_i[1]

Inputs are read from the datablock labels set by each C++ integrand's
module_label(). The 1-halo ΔΣ has an R axis supplied in the ini so we can
reshape the flat integration output back to (n_bins, n_r).
"""

import numpy as np
from cosmosis.datablock import option_section


def setup(options):
    try:
        nc_label = options.get_string(option_section, "num_counts_label")
    except Exception:
        nc_label = "NumCountsSel"
    try:
        mass_label = options.get_string(option_section, "mass_label")
    except Exception:
        mass_label = "MassWeightedSel"
    try:
        bias_label = options.get_string(option_section, "bias_label")
    except Exception:
        bias_label = "BiasWeightedSel"
    try:
        dsig_label = options.get_string(option_section, "dsigma_1h_label")
    except Exception:
        dsig_label = "DSigma1hSel"
    try:
        out_section = options.get_string(option_section, "output_section")
    except Exception:
        out_section = "cluster_observables"
    try:
        r_perp = np.asarray(
            options[option_section, "r_perp"], dtype=np.float64
        )
    except Exception:
        r_perp = None

    return {
        "nc_label": nc_label,
        "mass_label": mass_label,
        "bias_label": bias_label,
        "dsig_label": dsig_label,
        "out_section": out_section,
        "r_perp": r_perp,
    }


def _flat_vals(block, label):
    v = np.asarray(block[label, "vals"], dtype=np.float64)
    return v.reshape(-1)


def execute(block, config):
    out = config["out_section"]

    nc = _flat_vals(block, config["nc_label"])          # (n_bins,)
    mass = _flat_vals(block, config["mass_label"])      # (n_bins,)
    bias = _flat_vals(block, config["bias_label"])      # (n_bins,)
    dsig = _flat_vals(block, config["dsig_label"])      # (n_bins * n_r,)

    if nc.shape != mass.shape or nc.shape != bias.shape:
        raise ValueError(
            f"Shape mismatch in N-operator outputs: "
            f"nc={nc.shape}, mass={mass.shape}, bias={bias.shape}"
        )

    n_bins = nc.size
    m_eff = mass / nc
    b_eff = bias / nc

    block[out, "n_counts"] = nc
    block[out, "m_eff"] = m_eff
    block[out, "b_eff"] = b_eff

    if dsig.size % n_bins != 0:
        raise ValueError(
            f"dsigma_1h size {dsig.size} is not a multiple of n_bins {n_bins}"
        )
    n_r = dsig.size // n_bins
    dsig_2d = dsig.reshape(n_bins, n_r)
    dsig_1h = dsig_2d / nc[:, None]
    block[out, "delta_sigma_1h"] = dsig_1h
    if config["r_perp"] is not None:
        block[out, "r_perp"] = config["r_perp"]
    return 0


def cleanup(config):
    return 0
