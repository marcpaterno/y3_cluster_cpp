"""CosmoSIS Python module: replace CAMB's linear P(k) computation with
CosmoPower NN emulators loaded as torch-free / TF-free .npz artifacts.

Writes (subset of) `matter_power_lin`, `cdm_baryon_power_lin`, and
`matter_power_nl` sections to the datablock, matching the grid layout
that `mf_tinker` and `sigma_cpp` consume:

    block.put_grid(section, "z", z, "k_h", k, "p_k", p_k)
    p_k.shape == (n_z, n_k)

Also publishes `distances/z` and `distances/d_a` via astropy's
FlatLambdaCDM cosmology so downstream modules that depend on angular
diameter distances (e.g. DV_DO_DZ_t in y3_cluster_cpp) work without
CAMB. This can be turned off by setting `write_distances = false`.

ini options (in [cp_camb] section):
    emulator_repo         : path to the camb-emulator repo root so
                             `cp_numpy` is importable (required unless
                             already on PYTHONPATH).
    linear_pk_path        : .npz for delta_tot linear P(k) (REQUIRED).
                             Published to `matter_power_lin`.
    linear_nonu_pk_path   : .npz for delta_nonu linear P(k) (optional).
                             Published to `cdm_baryon_power_lin`.
    nonlinear_pk_path     : .npz for delta_tot non-linear P(k) (optional).
                             Published to `matter_power_nl`.
    nonu_fallback         : bool. If true and linear_nonu_pk_path is
                             not supplied, the linear emulator's output
                             is ALSO written to `cdm_baryon_power_lin`
                             so `mf_tinker matter_power_lin_version=2`
                             has something to read. Temporary hack.
    zmin, zmax, nz        : output z-grid (default 0.0, 1.0, 50).
    write_distances       : bool (default true). If true, also publish
                             `distances/{z,d_a,d_m,d_l,h,mu}` via
                             astropy FlatLambdaCDM.
    apply_growth          : bool (default true when an emulator has no
                             `z` parameter, false otherwise). When true,
                             reads `growth_parameters/{z, d_z}` from the
                             datablock and scales P(k, z=0) up to
                             P(k, z) = D(z)² · P(k, 0). Requires
                             `structure/growth_factor` to run BEFORE
                             cp_camb in the pipeline.
"""

from __future__ import annotations

import logging
import os
import sys

import numpy as np
from cosmosis.datablock import names, option_section

logger = logging.getLogger("cp_camb")
if not logger.handlers:
    logger.setLevel(logging.INFO)
    _h = logging.StreamHandler()
    _h.setFormatter(logging.Formatter("[cp_camb] %(levelname)s: %(message)s"))
    logger.addHandler(_h)

cosmo = names.cosmological_parameters

# CosmoPower emulators were trained on this parameter order, and
# parameter names are lowercased by CosmoSIS when reading the datablock.
_COSMO_PARAM_MAP = {
    "h0": "h0",
    "omega_m": "omega_m",
    "omega_b": "omega_b",
    "n_s": "n_s",
    "log1e10As": "log1e10as",   # CosmoSIS lowercases; be defensive
    "mnu": "mnu",
}


def _ensure_repo_on_path(options):
    try:
        repo = options.get_string(option_section, "emulator_repo")
    except Exception:
        repo = None
    if repo and repo not in sys.path:
        sys.path.insert(0, repo)


def _get_opt_path(options, key: str):
    try:
        p = options.get_string(option_section, key)
        return p if p else None
    except Exception:
        return None


def _get_opt_float(options, key: str, default: float) -> float:
    try:
        return float(options.get_double(option_section, key))
    except Exception:
        try:
            return float(options.get_int(option_section, key))
        except Exception:
            return default


def _get_opt_int(options, key: str, default: int) -> int:
    try:
        return int(options.get_int(option_section, key))
    except Exception:
        return default


def _get_opt_bool(options, key: str, default: bool) -> bool:
    try:
        return bool(options.get_bool(option_section, key))
    except Exception:
        return default


def setup(options):
    _ensure_repo_on_path(options)
    from cp_numpy import CosmoPowerNumpyNN

    linear_path = _get_opt_path(options, "linear_pk_path")
    if linear_path is None:
        raise ValueError(
            "cp_camb: linear_pk_path is required in the [cp_camb] section"
        )
    linear_nonu_path = _get_opt_path(options, "linear_nonu_pk_path")
    nonlinear_path = _get_opt_path(options, "nonlinear_pk_path")
    nonu_fallback = _get_opt_bool(options, "nonu_fallback", False)

    zmin = _get_opt_float(options, "zmin", 0.0)
    zmax = _get_opt_float(options, "zmax", 1.0)
    nz = _get_opt_int(options, "nz", 50)
    z_array = np.linspace(zmin, zmax, nz, dtype=np.float64)
    write_distances = _get_opt_bool(options, "write_distances", True)

    linear = CosmoPowerNumpyNN(linear_path)
    nonu = CosmoPowerNumpyNN(linear_nonu_path) if linear_nonu_path else None
    nonlin = CosmoPowerNumpyNN(nonlinear_path) if nonlinear_path else None

    # Detect whether any loaded emulator lacks `z` (v2-style, trained at
    # z=0 only). If so, we need D(z) to reconstruct P(k, z) = D(z)²·P(k,0).
    any_z0_only = any(
        emu is not None and "z" not in emu.parameters
        for emu in (linear, nonu, nonlin)
    )
    apply_growth = _get_opt_bool(options, "apply_growth", any_z0_only)
    if apply_growth and not any_z0_only:
        logger.warning(
            "apply_growth=True but all emulators have `z` as an input; "
            "growth scaling will be ignored."
        )
    if any_z0_only and not apply_growth and zmax > zmin:
        logger.warning(
            "At least one emulator is z=0-only but apply_growth=False. "
            "P(k, z>0) will be returned as P(k, 0) — not physically "
            "correct unless zmax=zmin=0. Chain "
            "`structure/growth_factor` before cp_camb and set "
            "apply_growth=True to fix."
        )

    # Sanity-check that all loaded emulators share the same k-modes; if
    # they don't, downstream consumers would see inconsistent grids.
    for label, emu in (("linear_nonu", nonu), ("nonlinear", nonlin)):
        if emu is None:
            continue
        if emu.modes.shape != linear.modes.shape or not np.allclose(
            emu.modes, linear.modes, rtol=1e-6
        ):
            raise ValueError(
                f"cp_camb: {label} emulator k-modes differ from linear; "
                "all emulators must share the same k-grid"
            )

    logger.info("loaded linear_pk_path=%s", os.path.basename(linear_path))
    if nonu is not None:
        logger.info("loaded linear_nonu_pk_path=%s",
                    os.path.basename(linear_nonu_path))
    elif nonu_fallback:
        logger.warning(
            "nonu_fallback=True: publishing linear emulator output to BOTH "
            "matter_power_lin AND cdm_baryon_power_lin. This is a "
            "placeholder until the real linear_nonu emulator is trained."
        )
    if nonlin is not None:
        logger.info("loaded nonlinear_pk_path=%s",
                    os.path.basename(nonlinear_path))
    logger.info("z grid: %.3f..%.3f (%d points)", zmin, zmax, nz)

    return {
        "linear": linear,
        "nonu": nonu,
        "nonlin": nonlin,
        "nonu_fallback": nonu_fallback,
        "z_array": z_array,
        "k_modes": linear.modes,
        "write_distances": write_distances,
        "apply_growth": apply_growth,
    }


def _read_cosmo_params(block):
    params = {}
    for emu_key, block_key in _COSMO_PARAM_MAP.items():
        # CosmoSIS lowercases by default; be tolerant of either case.
        if block.has_value(cosmo, block_key):
            params[emu_key] = float(block[cosmo, block_key])
        elif block.has_value(cosmo, emu_key.lower()):
            params[emu_key] = float(block[cosmo, emu_key.lower()])
        else:
            raise KeyError(
                f"cp_camb: {block_key!r} not found in cosmological_parameters"
            )
    return params


def _evaluate(emu, cosmo_params, z_array, growth=None):
    """Returns P(k, z) shape (n_z, n_k) in linear units (not log10).

    If `growth` is an array of shape (n_z,) matching z_array, evaluate
    the emulator once at z=0 and broadcast P(k, z) = D(z)² · P(k, 0).
    Use this path for z=0-only emulators (v2 CosmoPower: `z` not in
    the parameter list). Otherwise evaluate the emulator on the full
    z grid as before.
    """
    if growth is not None:
        log10_pk0 = emu.predictions_np_at_z(cosmo_params, np.array([0.0]))
        p0 = np.power(10.0, log10_pk0[0], dtype=np.float64)  # (n_k,)
        return (growth.astype(np.float64)[:, None] ** 2) * p0[None, :]
    log10_pk = emu.predictions_np_at_z(cosmo_params, z_array)
    return np.power(10.0, log10_pk, dtype=np.float64)


def _read_growth(block, z_array):
    """Return D(z) on `z_array` interpolated from
    `growth_parameters/{z, d_z}`, renormalized so that D(z=0) = 1.

    The CosmoSIS `structure/growth_factor` module normalizes D to the
    matter-dominated asymptote (D(z_early) = a), so D(0) ≈ 0.76 for
    ΛCDM — NOT 1. The camb-emulator was trained on P(k, z=0), so to
    reconstruct P(k, z) = D(z)² · P(k, 0) consistently we divide out
    the raw D(0).

    Raises if the section is missing — the caller is expected to chain
    `structure/growth_factor` ahead of cp_camb.
    """
    from cosmosis.datablock import names as _names
    section = _names.growth_parameters
    if not block.has_value(section, "z") or not block.has_value(section, "d_z"):
        raise KeyError(
            "cp_camb: apply_growth=True but growth_parameters/{z,d_z} not "
            "on the datablock. Add `structure/growth_factor` to the "
            "pipeline before cp_camb."
        )
    z_in = np.asarray(block[section, "z"], dtype=np.float64)
    d_in = np.asarray(block[section, "d_z"], dtype=np.float64)
    d_at_0 = float(np.interp(0.0, z_in, d_in))
    if d_at_0 <= 0:
        raise ValueError(
            f"cp_camb: invalid D(z=0)={d_at_0} from growth_parameters"
        )
    return np.interp(z_array, z_in, d_in) / d_at_0


def _write_distances(block, cosmo_params, z_array):
    """Emit `distances/{z, d_a, d_m, d_l, h, mu}` via astropy's
    FlatLambdaCDM. Matches the convention used by
    y3_buzzard/halo_model_cosmosis.py:59 and the default pipeline values
    expected by dv_do_dz_t.hh (distances in h⁻¹ Mpc — but CAMB outputs
    distances in Mpc; we match CAMB's convention).
    """
    import astropy.cosmology

    h0 = cosmo_params["h0"]
    H0_kms = 100.0 * h0
    omega_m = cosmo_params["omega_m"]
    omega_b = cosmo_params["omega_b"]
    cosmo = astropy.cosmology.FlatLambdaCDM(
        H0=H0_kms,
        Om0=omega_m,
        Ob0=omega_b,
        Tcmb0=2.725,
    )

    # astropy returns astropy Quantities; drop units to plain floats.
    z = np.asarray(z_array, dtype=np.float64)
    # Avoid z=0 at the front because d_a→0 there; shift epsilon if needed.
    # d_a in Mpc (NOT h⁻¹): CAMB's convention.
    d_a = cosmo.angular_diameter_distance(z).to_value("Mpc")
    d_m = cosmo.comoving_transverse_distance(z).to_value("Mpc")
    d_l = cosmo.luminosity_distance(z).to_value("Mpc")
    H_z = cosmo.H(z).to_value("km / (Mpc s)")
    # mu: distance modulus (CAMB units: magnitudes). Guard z=0.
    with np.errstate(divide="ignore"):
        mu = 5.0 * np.log10(d_l) + 25.0
        mu = np.where(np.isfinite(mu), mu, 0.0)

    block["distances", "z"] = z
    block["distances", "a"] = 1.0 / (1.0 + z)
    block["distances", "d_a"] = d_a
    block["distances", "d_m"] = d_m
    block["distances", "d_l"] = d_l
    # Flat-cosmology comoving radial distance equals d_m; emit under the
    # `d_c` name that y3_cluster_cpp's P_operator / sigma_prj read.
    block["distances", "d_c"] = d_m
    # CAMB convention: distances/h = H(z)/c in 1/Mpc
    # (c = 299792.458 km/s; H_z in km/s/Mpc).
    block["distances", "h"] = H_z / 299792.458
    block["distances", "mu"] = mu
    block["distances", "nz"] = len(z)


def execute(block, config):
    cosmo_params = _read_cosmo_params(block)
    z = config["z_array"]
    k = config["k_modes"]

    if config["write_distances"]:
        _write_distances(block, cosmo_params, z)

    # D(z) once per sample; reused by every z=0-only emulator below.
    d_z = _read_growth(block, z) if config["apply_growth"] else None

    def _growth_for(emu):
        # Only apply growth scaling to emulators that lack `z` as an input.
        return d_z if (config["apply_growth"] and "z" not in emu.parameters) else None

    p_lin = _evaluate(config["linear"], cosmo_params, z, _growth_for(config["linear"]))
    block.put_grid("matter_power_lin", "z", z, "k_h", k, "p_k", p_lin)

    if config["nonu"] is not None:
        p_nonu = _evaluate(
            config["nonu"], cosmo_params, z, _growth_for(config["nonu"])
        )
        block.put_grid(
            "cdm_baryon_power_lin", "z", z, "k_h", k, "p_k", p_nonu
        )
    elif config["nonu_fallback"]:
        # Publish the linear total-matter P(k) under the nonu section too,
        # so mf_tinker matter_power_lin_version=2 has a grid to read.
        block.put_grid(
            "cdm_baryon_power_lin", "z", z, "k_h", k, "p_k", p_lin
        )

    if config["nonlin"] is not None:
        p_nl = _evaluate(
            config["nonlin"], cosmo_params, z, _growth_for(config["nonlin"])
        )
        block.put_grid("matter_power_nl", "z", z, "k_h", k, "p_k", p_nl)

    return 0


def cleanup(config):
    return 0
