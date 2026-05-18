"""Combine Sigma_prj, DSigma_prj, and Sigma_crit^-1 into reduced shear
for the projection-effects two-halo term.

Inputs (datablock):
    sigma_prj/vals       (n_lob * n_zob * n_R,)   <Sigma_prj>
    dsigma_prj/vals      (n_lob * n_zob * n_R,)   <DSigma_prj>
    average_sigma_crit_inv/{zlense, sci_average}  Sigma_crit^-1(z_lens)

Output (datablock):
    red_shear_prj/vals   (n_lob * n_zob * n_R,)   g_t^prj
    red_shear_prj/gamma_t_prj                     gamma_t (numerator only)

Formula (per grid point (lob_i, zob_j, R_k)):
    gamma_t_prj = <DSigma_prj> * Sigma_crit^-1(z_lens = (zob_lo + zob_hi)/2)
    g_t^prj     = gamma_t_prj / (1 - <Sigma_prj> * Sigma_crit^-1)

Grid axes come from this module's own ini section:
    zo_low, zo_high, lambda_bin, radii  (all length n_lob * n_zob * n_R,
                                         zipped as (zob_idx, lambda_idx, R_idx)
                                         in the same order as sigma_prj's
                                         wall-of-numbers grid).
"""
from __future__ import annotations
import numpy as np
from cosmosis.datablock import option_section


def _doubles(options, key, required=True):
    try:
        return np.asarray(
            options.get_double_array_1d(option_section, key), dtype=float)
    except Exception:
        try:
            return np.asarray(
                options.get_int_array_1d(option_section, key), dtype=float)
        except Exception:
            if required:
                raise
            return None


def setup(options):
    cfg = {}
    cfg['zo_low']     = _doubles(options, 'zo_low')
    cfg['zo_high']    = _doubles(options, 'zo_high')
    cfg['lambda_bin'] = _doubles(options, 'lambda_bin')
    cfg['radii']      = _doubles(options, 'radii')

    # Optional: input/output section names so this can be reused for
    # other (Sigma, DSigma) pairs.
    try:
        cfg['sigma_section'] = options.get_string(option_section,
                                                   'sigma_section')
    except Exception:
        cfg['sigma_section'] = 'sigma_prj'
    try:
        cfg['dsigma_section'] = options.get_string(option_section,
                                                    'dsigma_section')
    except Exception:
        cfg['dsigma_section'] = 'dsigma_prj'
    try:
        cfg['output_section'] = options.get_string(option_section,
                                                    'output_section')
    except Exception:
        cfg['output_section'] = 'red_shear_prj'
    try:
        cfg['sci_section'] = options.get_string(option_section,
                                                 'sci_section')
    except Exception:
        cfg['sci_section'] = 'average_sigma_crit_inv'

    n = cfg['zo_low'].size
    for k in ('zo_high', 'lambda_bin', 'radii'):
        if cfg[k].size != n:
            raise ValueError(
                f"red_shear_prj: axis '{k}' length {cfg[k].size} != {n}")
    cfg['n_grid'] = n
    return cfg


def execute(block, config):
    import time
    t0 = time.perf_counter()

    # Flat arrays, same order as the Sigma_prj / DSigma_prj evaluators wrote.
    sig  = np.asarray(block[config['sigma_section'], 'vals'], dtype=float).ravel()
    dsig = np.asarray(block[config['dsigma_section'], 'vals'], dtype=float).ravel()
    if sig.shape != dsig.shape or sig.size != config['n_grid']:
        raise ValueError(
            f"red_shear_prj: vals shape mismatch "
            f"sigma={sig.shape} dsigma={dsig.shape} expected {config['n_grid']}")

    # Sigma_crit^-1(z_lens) at the midpoint of each (zo_low, zo_high) bin.
    z_sci = np.asarray(block[config['sci_section'], 'zlense'], dtype=float)
    sci_v = np.asarray(block[config['sci_section'], 'sci_average'], dtype=float)
    z_mid = 0.5 * (config['zo_low'] + config['zo_high'])
    sci_at_z = np.interp(z_mid, z_sci, sci_v)            # (n_grid,)

    gamma_t = dsig * sci_at_z
    denom   = 1.0 - sig * sci_at_z
    # Guard against pathological (1 - Σ·Σ_cr^-1) ≤ 0 (happens inside
    # critical-density source distributions; not a physical issue for
    # cluster cosmology but protect the ratio numerically).
    g_t = np.where(denom > 1.0e-12, gamma_t / denom, gamma_t)

    out = config['output_section']
    block[out, 'vals']        = g_t
    block[out, 'gamma_t_prj'] = gamma_t
    block[out, 'zo_low']      = config['zo_low']
    block[out, 'zo_high']     = config['zo_high']
    block[out, 'lambda_bin']  = config['lambda_bin']
    block[out, 'radii']       = config['radii']

    dt_ms = 1000.0 * (time.perf_counter() - t0)
    print(f"[red_shear_prj] assembled g_t^prj on {config['n_grid']} "
          f"grid points in {dt_ms:.1f} ms", flush=True)
    return 0


def cleanup(config):
    return 0
