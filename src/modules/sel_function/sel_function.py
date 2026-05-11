"""Tabulate the richness-selection function S_ij(lnM, z) once per sample.

Publishes a single packed tensor on a shared (lnM, z) grid that all 12
wall-grid bins live on. The downstream C++ integrand (SelFunction_t)
slices the plane for its bin and serves it via Interp2D.

    S_ij(lnM, z) = S_i(lnM, z) * K_j(z)
                 = [ sum_k W_k K_i(lam_k, z) P_HOD(lam_k | M, z) ] * K_j(z)

with
    K_i(ltr, z)  — closed-form bin-integrated richness kernel (plob_ltr EMG).
                   Computed on the fly by CDF differencing at unique bin
                   edges (e.g. {20, 30, 45, 60, 200}) — no table, no cache.
    K_j(z)       — Gaussian CDF difference over the z bin.
    P_HOD        — Shifted-Poisson continuous form (doc Eq. 28), matching
                   src/models/mor_hod_t.hh. Single gammaln call on the
                   (n_lnm, n_z, N_q) tensor.

Grid:
    lnM  — linear in [lnm_low, lnm_high], n_lnm points (default 256).
    z    — linear in [zt_low, zt_high], n_z_shared points (default 64).
           Individual bins hard-zero outside |ztr - zob_mid| > L_z sigma_z
           via the K_j factor.

Reads (ini options):
    sel_function/{lam_min, lam_max, zob_min, zob_max, sigma_z}   per-bin
    sel_function/{zt_low, zt_high, lnm_low, lnm_high}            envelope
    sel_function/{n_lnm, n_z_shared, L_z, N_q}                   grid knobs

Reads (datablock):
    cluster_mor/{log10_Mmin, log10_M1|log10_ratio, alpha, epsilon,
                 sigma_lambda, z_pivot?}
    plob_ltr_params/{z, a_tau, b_tau, a_mu, b_mu, a_sig, b_sig,
                     a_fprj, b_fprj}

Writes to datablock (section "sel_function"):
    lnM       (n_lnm,)             shared mass grid
    z         (n_z,)               shared redshift grid
    S_stack   (n_bins, n_z, n_lnm) S_ij for every bin in one tensor
"""
from __future__ import annotations
import os
import numpy as np
from cosmosis.datablock import option_section
from scipy.special import erf, erfcx, gammaln


SQRT2 = np.sqrt(2.0)
SQRT2PI = np.sqrt(2.0 * np.pi)


# ---------------------------------------------------------------------------
# Ini reading helpers
# ---------------------------------------------------------------------------

def _read_array(options, section, name, dtype=float):
    try:
        return np.asarray(options.get_double_array_1d(section, name), dtype=dtype)
    except Exception:
        return np.asarray(options.get_int_array_1d(section, name), dtype=dtype)


def _read_scalar_or_first(options, section, name, default):
    try:
        v = options.get_double_array_1d(section, name)
        return float(np.asarray(v, dtype=float).ravel()[0])
    except Exception:
        try:
            return float(options.get_double(section, name))
        except Exception:
            return float(default)


# ---------------------------------------------------------------------------
# Closed-form pieces
# ---------------------------------------------------------------------------

def _phi(x):
    """Standard normal CDF."""
    return 0.5 * (1.0 + erf(x / SQRT2))


def _f_emg(x, mu, sigma, tau):
    """EMG CDF via erfcx — broadcasts over ndarrays.

    Matches src/models/richness_kernel_t.hh::F_EMG.
    """
    z = (x - mu) / sigma
    u = (tau * sigma - z) / SQRT2
    neg = u < 0.0
    abs_u = np.where(neg, -u, u)
    exp_mz2 = np.exp(-0.5 * z * z)
    tail_base = 0.5 * erfcx(abs_u) * exp_mz2
    A = -tau * (x - mu) + 0.5 * (tau * sigma) ** 2
    tail = np.where(neg, np.exp(np.clip(A, -700.0, 700.0)) - tail_base,
                    tail_base)
    return np.clip(_phi(z) - tail, 0.0, 1.0)


def _plob_params(ltr, z, plob_splines):
    """Evaluate the 4 EMG params (mu, sigma, tau, fprj) at (ltr, z).

    Centralises the (ltr, z) closed-form costs so the per-bin K_i loop
    can reuse the output via CDF differencing (see _K_edges_of_bins).

    Fast path: when ``z`` is a 1-D array of shape ``(N_z,)`` and ``ltr``
    has shape ``(N_lnm, N_z, N_q)`` (the production call site), each of
    the 8 coefficient splines is evaluated on the 1-D grid only (64
    points) instead of the full 524k-element broadcast tensor — saves
    ~130 ms/sample at ``(256, 64, 32)``.  Fallback path (same ndim as
    ltr) is kept for the diagnostic helpers.
    """
    z_arr = np.asarray(z, dtype=float)
    ltr_ndim = np.ndim(ltr)
    if z_arr.ndim == 1 and ltr_ndim > 1:
        # Broadcast axis: for (n_lnm, n_z, n_q) this is axis=-2 (the z axis);
        # reshape the per-z coefficients to (1, N_z, 1, ..., 1).  Other
        # layouts fall back through a more generic reshape search.
        target_axis = None
        shape_ltr = np.shape(ltr)
        for ax, s in enumerate(shape_ltr):
            if s == z_arr.size:
                target_axis = ax
                break
        if target_axis is None:
            # Give up on the fast path — shapes don't agree on which axis is z.
            bshape = ()
        else:
            bshape = tuple(z_arr.size if ax == target_axis else 1
                           for ax in range(ltr_ndim))
        def _eval(k):
            v = plob_splines[k](z_arr)
            return v.reshape(bshape) if bshape else v
    else:
        def _eval(k):
            return plob_splines[k](z_arr)

    a_mu   = _eval('a_mu')
    b_mu   = _eval('b_mu')
    a_sig  = _eval('a_sig')
    b_sig  = _eval('b_sig')
    a_tau  = _eval('a_tau')
    b_tau  = _eval('b_tau')
    a_fprj = _eval('a_fprj')
    b_fprj = _eval('b_fprj')

    mu    = a_mu + b_mu * ltr
    sigma = b_sig * np.power(ltr, a_sig)
    tau   = b_tau / np.power(ltr, a_tau)
    fprj  = np.minimum(1.0, b_fprj / np.power(1.0 + np.exp(-ltr), a_fprj))
    return mu, sigma, tau, fprj


def _cdf_lob(lam_ob, mu, sigma, tau, fprj):
    """CDF of the Costanzi EMG kernel at lam_ob, given (mu, sigma, tau, fprj)
    already evaluated at (ltr, z).

    CDF(lob | ltr, z) = Phi(z_std) - fprj * tail
        tail = exp(-tau (lob-mu) + tau^2 sigma^2 / 2) * Phi(tau sigma - z_std)
             = 0.5 * erfcx(|u|) * exp(-0.5 z_std^2)            (u >= 0)
             = exp(A) - 0.5 * erfcx(|u|) * exp(-0.5 z_std^2)    (u <  0)
    with u = (tau sigma - z_std)/sqrt(2),  A = -tau(lob-mu) + 0.5(tau sigma)^2.

    `lam_ob` is a scalar; `mu, sigma, tau, fprj` broadcast (same shape).

    Optimisation notes: scalar kernel path via `np.subtract(out=)` and
    in-place updates to the single pre-allocated `tail` buffer.  The
    critical change vs the previous version is dropping the branched
    `np.where(neg, exp(clip(A)) - tail_base, tail_base)` which
    evaluated `np.exp(np.clip(A, ...))` on all 524k elements even where
    the `u>=0` branch was taken; we now only compute `exp(A)` for the
    `neg` subset.
    """
    z_std = (lam_ob - mu) / sigma
    u = (tau * sigma - z_std) * (1.0 / SQRT2)
    # tail_base = 0.5 * erfcx(|u|) * exp(-0.5 z_std^2).  abs(u) instead of
    # np.where(neg, -u, u) -- faster and same result since erfcx is even
    # for our sign convention below.
    exp_mz2 = np.exp(-0.5 * z_std * z_std)
    tail = 0.5 * erfcx(np.abs(u)) * exp_mz2                    # shape of mu
    # Correction in the u<0 branch: tail -> exp(A) - tail_base.
    # Compute exp(A) ONLY on the negative-u subset (saves 524k exps).
    neg = u < 0.0
    if np.any(neg):
        A_neg = -tau[neg] * (lam_ob - mu[neg]) + 0.5 * (tau[neg] * sigma[neg]) ** 2
        np.clip(A_neg, -700.0, 700.0, out=A_neg)
        # tail[neg] = exp(A_neg) - tail[neg]
        tail[neg] = np.exp(A_neg) - tail[neg]
    # CDF_Costanzi = Phi(z_std) - fprj * tail.  fuse subtract + clip.
    out = _phi(z_std)
    out -= fprj * tail
    np.clip(out, 0.0, 1.0, out=out)
    return out


def _cdf_lob_stacked(lam_edges, mu, sigma, tau, fprj):
    """Vectorised form of _cdf_lob over a 1-D set of lam_ob edges.

    Pre-computes ``tau*sigma``, ``(tau*sigma)^2/2``, ``1/sigma`` once
    across all edges (small gain, ~5 ms) and keeps the per-edge hot
    path identical to _cdf_lob.  A true ``(n_edges, *shape)`` stacked
    ufunc call was tried and found to be SLOWER than the per-edge
    loop — scipy's erf/erfcx are memory-bandwidth-limited and stacking
    5 copies into one call costs more in allocation than it saves in
    dispatch overhead.  Keep the loop + shared arithmetic only.

    Matches ``_cdf_lob`` bit-for-bit at every edge.
    """
    tau_sigma = tau * sigma
    half_ts2  = 0.5 * tau_sigma * tau_sigma
    inv_sigma = 1.0 / sigma
    outs = []
    for lam_ob in lam_edges:
        lam = float(lam_ob)
        z_std = (lam - mu) * inv_sigma
        u = (tau_sigma - z_std) * (1.0 / SQRT2)
        tail = 0.5 * erfcx(np.abs(u)) * np.exp(-0.5 * z_std * z_std)
        neg = u < 0.0
        if np.any(neg):
            A_neg = -tau[neg] * (lam - mu[neg]) + half_ts2[neg]
            np.clip(A_neg, -700.0, 700.0, out=A_neg)
            tail[neg] = np.exp(A_neg) - tail[neg]
        out = _phi(z_std)
        out -= fprj * tail
        np.clip(out, 0.0, 1.0, out=out)
        outs.append(out)
    return outs


def _K_edges_of_bins(lam_edges, ltr, z, plob_splines):
    """Compute K_i for every bin in one shot by differencing the CDF at the
    shared bin edges.

    lam_edges : 1-D array of unique bin edges in increasing order, length
                n_bins + 1. Example Y3: [20, 30, 45, 60, 200].

    Returns (cdf, K_per_bin) where
        cdf       has shape (*ltr.shape, n_edges)
        K_per_bin has shape (*ltr.shape, n_bins)
    """
    mu, sigma, tau, fprj = _plob_params(ltr, z, plob_splines)
    # Stack CDF over edges — each edge share mu/sigma/tau/fprj.
    cdfs = np.stack(
        [_cdf_lob(x_e, mu, sigma, tau, fprj) for x_e in lam_edges],
        axis=-1,
    )
    # K_i = CDF(edge_{i+1}) - CDF(edge_i)
    K_per_bin = cdfs[..., 1:] - cdfs[..., :-1]
    return cdfs, K_per_bin


def _K_i_bin(lam_min, lam_max, ltr, z, plob_splines):
    """K_i(ltr, z) closed form — kept for parity / debug. Production path
    goes through _K_edges_of_bins for the full bin set."""
    mu, sigma, tau, fprj = _plob_params(ltr, z, plob_splines)
    cdf_hi = _cdf_lob(lam_max, mu, sigma, tau, fprj)
    cdf_lo = _cdf_lob(lam_min, mu, sigma, tau, fprj)
    return cdf_hi - cdf_lo


def _K_j(ztr, zob_min, zob_max, sigma_z):
    """Gaussian CDF difference over the z-bin."""
    return _phi((zob_max - ztr) / sigma_z) - _phi((zob_min - ztr) / sigma_z)


# ---------------------------------------------------------------------------
# HOD MOR — Poisson(mu_sat) * Gauss(lambda_sigma * mu_sat), truncated at 0
# and renormalized. Matches src/models/mor_hod_t.hh line-for-line.
# ---------------------------------------------------------------------------

POISSON_TOL = 1e-8
FALLBACK_SIGMA = 1.0e-3
MIN_SIGMA = 1.0e-6
Z_PIVOT_DEFAULT = 0.45


def _mu_sat(M, z, log10_Mmin, log10_M1, alpha, epsilon, z_pivot):
    """HOD mean satellite count as a function of (M, z). Broadcasts."""
    Mmin = 10.0 ** log10_Mmin
    M1   = 10.0 ** log10_M1
    dM1  = M1 - Mmin
    if dM1 <= 0.0:
        return np.zeros_like(M)
    dM = np.maximum(M - Mmin, 0.0)
    base = np.where(dM > 0.0, dM / dM1, 0.0)
    with np.errstate(invalid='ignore', divide='ignore'):
        mu = np.where(
            base > 0.0,
            np.power(base, alpha) * np.power((1.0 + z) / (1.0 + z_pivot),
                                             epsilon),
            0.0,
        )
    return mu


def _p_hod_scalar(ltr, lnM, z, mor):
    """P_HOD(ltr | M, z) — shifted continuous-Poisson form (doc Eq. 28).

    With nu = mu_sat + delta and delta = (sigma_intr * mu_sat)^2:

        P(ltr | M, z) = exp[ -nu + (ltr + delta - 1) ln(nu) - lnGamma(ltr + delta) ]

    Closed-form, no discrete k-sum, no 4-D (k) tensor. Numpy vectorised
    across the full (n_lnm, n_z, N_q) ltr tensor in a single gammaln call.

    ltr : (..., N_q)    — lambda quadrature nodes for each (M, z) cell
    lnM : (n_lnm, 1, 1) broadcasts along z, ltr
    z   : (1, n_z, 1)
    mor : dict of scalar HOD params
    """
    M = np.exp(lnM)                                            # (n_lnm, 1, 1)
    mu_sat = _mu_sat(M, z,
                     mor['log10_Mmin'], mor['log10_M1'],
                     mor['alpha'], mor['epsilon'], mor['z_pivot'])
    lcentral = np.where(M >= 10.0 ** mor['log10_Mmin'], 1.0, 0.0)

    # delta = (sigma_intr * mu_sat)^2;  nu = mu_sat + delta  (both (n_lnm, n_z, 1))
    delta = (mor['sigma_lambda'] * mu_sat) ** 2
    nu = mu_sat + delta

    # Where mu_sat is tiny, the HOD reduces to a delta at lambda_central
    # (represented as a narrow Gaussian for integration stability).
    tiny = mu_sat <= POISSON_TOL
    fallback = (np.exp(-0.5 * ((ltr - lcentral) / FALLBACK_SIGMA) ** 2)
                / (SQRT2PI * FALLBACK_SIGMA))

    # Continuous shifted Poisson. Shift argument by (lambda_central + delta)
    # so the "mean" sits at lcentral + mu_sat (matches doc Eq. 28 where
    # ltr is the latent richness including the central).
    x = ltr - lcentral + delta                              # broadcast to full tensor
    log_nu = np.log(np.maximum(nu, 1e-300))
    # Only valid for x > 0 (i.e. ltr > lcentral - delta); clip and zero out.
    valid = x > 0.0
    x_safe = np.where(valid, x, 1.0)
    log_P = -nu + (x_safe - 1.0) * log_nu - gammaln(x_safe)
    density = np.where(valid, np.exp(log_P), 0.0)

    # Collapse to narrow Gaussian where mu_sat ~ 0 (limit of shifted Poisson
    # is ill-defined; matches the C++ fallback branch).
    density = np.where(tiny, fallback, density)
    # Physical boundaries: ltr >= 0 AND lambda_sat = ltr - lambda_central >= 0.
    # Zero out any ltr below lambda_central (or below 0).
    density = np.where((ltr < 0.0) | (ltr < lcentral), 0.0, density)
    return density


# ---------------------------------------------------------------------------
# Per-bin grid choice
# ---------------------------------------------------------------------------

def _choose_z_grid(zob_min, zob_max, sigma_z, zt_low, zt_high, L_z, n_z):
    """Tight z-grid where K_j has support, clipped to the integration volume."""
    z_lo = max(zt_low,  zob_min - L_z * sigma_z)
    z_hi = min(zt_high, zob_max + L_z * sigma_z)
    if z_hi <= z_lo:
        # Degenerate — tiny fallback grid (K_j is effectively zero everywhere)
        z_lo, z_hi = zt_low, zt_high
    return np.linspace(z_lo, z_hi, n_z)


def _choose_lnM_grid(lam_min, lam_max, zob_min, zob_max, mor,
                     lnm_low, lnm_high, n_lnm, lam_n_sigma=6.0):
    """Solve mu_sat(M, z_mid) for the bin's richness support.

    HOD scatter width at mean mu is sqrt(mu + (sigma_lambda*mu)**2).
    Lower: mu_sat ~ lam_min - L*sqrt(...); upper: mu_sat ~ lam_max + L*sqrt(...).
    Simpler robust choice: mu_sat in [lam_min/4, 2*lam_max] (covers ~6-sigma
    on the typical sigma_lambda ~ 0.2).
    """
    log10_Mmin = mor['log10_Mmin']
    log10_M1   = mor['log10_M1']
    alpha      = mor['alpha']
    epsilon    = mor['epsilon']
    z_pivot    = mor['z_pivot']
    z_mid = 0.5 * (zob_min + zob_max)

    Mmin = 10.0 ** log10_Mmin
    M1   = 10.0 ** log10_M1
    dM1  = M1 - Mmin
    redshift_evo = ((1.0 + z_mid) / (1.0 + z_pivot)) ** epsilon

    def M_of_mu(mu_target):
        # mu = ((M-Mmin)/dM1)^alpha * redshift_evo  =>
        # M = Mmin + dM1 * (mu/redshift_evo)^(1/alpha)
        return Mmin + dM1 * (max(mu_target / redshift_evo, 0.0)) ** (1.0 / alpha)

    # Widen vs GL-style [λ/4, 2·λ_max]: we need the whole HMF × S_i tail
    # to be captured, not just where S_i peaks. Empirically [λ/8, 4·λ_max]
    # holds bin-integrated error below ~0.1% vs the direct pipeline.
    mu_lo = max(0.125 * lam_min, 1e-3)
    mu_hi = 4.0 * lam_max if np.isfinite(lam_max) else None

    lnm_lo = np.log(max(M_of_mu(mu_lo), np.exp(lnm_low)))
    if mu_hi is None:
        lnm_hi = lnm_high
    else:
        lnm_hi = np.log(min(M_of_mu(mu_hi), np.exp(lnm_high)))
    if lnm_hi <= lnm_lo:
        lnm_lo, lnm_hi = lnm_low, lnm_high
    lnm_lo = max(lnm_lo, lnm_low)
    lnm_hi = min(lnm_hi, lnm_high)
    return np.linspace(lnm_lo, lnm_hi, n_lnm)


# ---------------------------------------------------------------------------
# Gauss–Legendre nodes (shared across bins, but bracket [a,b] is per-cell)
# ---------------------------------------------------------------------------

def _gl_nodes(N_q):
    t, w = np.polynomial.legendre.leggauss(N_q)
    return t, w


# ---------------------------------------------------------------------------
# CosmoSIS module entry points
# ---------------------------------------------------------------------------

def setup(options):
    cfg = {}
    cfg['lam_min'] = _read_array(options, option_section, 'lam_min')
    cfg['lam_max'] = _read_array(options, option_section, 'lam_max')
    cfg['zob_min'] = _read_array(options, option_section, 'zob_min')
    cfg['zob_max'] = _read_array(options, option_section, 'zob_max')
    cfg['sigma_z'] = _read_array(options, option_section, 'sigma_z')

    # Envelope of the downstream integration volume. Pass-through so the
    # feeder clips its per-bin grids to the same box.
    cfg['zt_low']   = _read_scalar_or_first(options, option_section, 'zt_low',
                                            0.05)
    cfg['zt_high']  = _read_scalar_or_first(options, option_section, 'zt_high',
                                            0.80)
    cfg['lnm_low']  = _read_scalar_or_first(options, option_section, 'lnm_low',
                                            np.log(1.0e13))
    cfg['lnm_high'] = _read_scalar_or_first(options, option_section, 'lnm_high',
                                            np.log(9.0e15))

    try:
        cfg['n_lnm'] = int(options.get_int(option_section, 'n_lnm'))
    except Exception:
        cfg['n_lnm'] = 40
    try:
        cfg['n_z'] = int(options.get_int(option_section, 'n_z'))
    except Exception:
        cfg['n_z'] = 20
    # Shared (lnM, z) grid size — all 12 bins live on this single grid.
    try:
        cfg['n_z_shared'] = int(options.get_int(option_section, 'n_z_shared'))
    except Exception:
        cfg['n_z_shared'] = cfg['n_z']
    try:
        cfg['L_z'] = float(options.get_double(option_section, 'L_z'))
    except Exception:
        cfg['L_z'] = 6.0
    try:
        cfg['N_q'] = int(options.get_int(option_section, 'N_q'))
    except Exception:
        cfg['N_q'] = 32
    try:
        cfg['L_lam'] = float(options.get_double(option_section, 'L_lam'))
    except Exception:
        cfg['L_lam'] = 6.0

    t, w = _gl_nodes(cfg['N_q'])
    cfg['gl_t'] = t
    cfg['gl_w'] = w

    n = cfg['lam_min'].size
    for key in ('lam_max', 'zob_min', 'zob_max', 'sigma_z'):
        if cfg[key].size != n:
            raise ValueError(
                f"sel_function: axis '{key}' size {cfg[key].size} != {n}")
    cfg['n_bins'] = n
    return cfg


# --- Plob spline helpers (linear in z, flat extrapolation past edges) ------

_PLOB_NPZ_DEFAULT = os.path.join(
    os.environ.get("Y3_CLUSTER_CPP_DIR",
                   "/pscratch/sd/j/jesteves/y3_cluster_cpp"),
    "data", "prj_params", "plob_ltr_params.npz",
)


def _make_plob_splines(block):
    """Return a dict of splines over z for the 8 EMG coefficients.

    Prefers the datablock section ``plob_ltr_params/*`` (populated
    by the legacy plob_ltr_loader module); falls back to the
    pre-baked npz at ``data/prj_params/plob_ltr_params.npz`` (the
    source of truth after plob_ltr_loader was retired 2026-05-06).
    """
    from scipy.interpolate import InterpolatedUnivariateSpline as ius
    keys = ('a_tau', 'b_tau', 'a_mu', 'b_mu', 'a_sig', 'b_sig',
            'a_fprj', 'b_fprj')
    try:
        z = np.asarray(block['plob_ltr_params', 'z'], dtype=float)
        src = {k: np.asarray(block['plob_ltr_params', k], dtype=float)
               for k in keys}
    except Exception:
        path = os.environ.get("RICHNESS_SELECTION_PLOB_NPZ",
                              _PLOB_NPZ_DEFAULT)
        data = np.load(path)
        z = np.asarray(data['z'], dtype=float)
        src = {k: np.asarray(data[k], dtype=float) for k in keys}
    return {k: ius(z, src[k], k=1, ext=3) for k in keys}


def _read_mor(block):
    """Read HOD parameters, supporting both log10_M1 and log10_ratio forms."""
    log10_Mmin = float(block['cluster_mor', 'log10_Mmin'])
    try:
        log10_ratio = float(block['cluster_mor', 'log10_ratio'])
        log10_M1 = log10_Mmin + log10_ratio
    except Exception:
        log10_M1 = float(block['cluster_mor', 'log10_M1'])
    z_pivot = Z_PIVOT_DEFAULT
    try:
        z_pivot = float(block['cluster_mor', 'z_pivot'])
    except Exception:
        pass
    return dict(
        log10_Mmin   = log10_Mmin,
        log10_M1     = log10_M1,
        alpha        = float(block['cluster_mor', 'alpha']),
        epsilon      = float(block['cluster_mor', 'epsilon']),
        sigma_lambda = float(block['cluster_mor', 'sigma_lambda']),
        z_pivot      = z_pivot,
    )


def _compute_lam_nodes_and_P_HOD(lnM, z, mor, gl_t, gl_w, L=6.0):
    """Per-(lnM, z) GL bracket + ltr nodes + P_HOD, bin-independent.

    Returns:
        lam_k       (n_lnm, n_z, n_q)   GL nodes in ltr
        W_k         (n_lnm, n_z, n_q)   GL weights scaled by (b-a)/2
        P_Mz        (n_lnm, n_z, n_q)   HOD P(ltr|M,z)
        degenerate  (n_lnm, n_z)        cells where b <= a (mu_sat ~ 0)
    """
    M = np.exp(lnM)[:, None]
    zz = z[None, :]
    mu_sat = _mu_sat(M, zz,
                     mor['log10_Mmin'], mor['log10_M1'],
                     mor['alpha'], mor['epsilon'], mor['z_pivot'])
    lcentral = np.where(M >= 10.0 ** mor['log10_Mmin'], 1.0, 0.0)
    mu_eff = lcentral + mu_sat
    sig_eff = np.sqrt(np.maximum(
        mu_sat + (mor['sigma_lambda'] * mu_sat) ** 2, 0.0))

    a = np.maximum(0.0, mu_eff - L * sig_eff)
    b = mu_eff + L * sig_eff
    degenerate = b <= a
    half = 0.5 * (b - a)
    mid  = 0.5 * (a + b)

    lam_k = mid[..., None] + half[..., None] * gl_t[None, None, :]
    W_k   = half[..., None] * gl_w[None, None, :]

    lnM_b = lnM[:, None, None]
    z_b   = z[None, :, None]
    P_Mz  = _p_hod_scalar(lam_k, lnM_b, z_b, mor)
    return lam_k, W_k, P_Mz, degenerate


def _unique_edges(lam_min, lam_max):
    """Sorted unique bin edges from per-bin (lam_min, lam_max)."""
    return np.unique(np.concatenate([lam_min, lam_max]))


def execute(block, config):
    import time
    t_start = time.perf_counter()

    mor = _read_mor(block)
    plob = _make_plob_splines(block)

    n_bins = config['n_bins']

    # Single shared (lnM, z) grid across all bins → compute P_HOD and the
    # GL bracket ONCE. Individual bins hard-zero outside their z-window in
    # K_j. The K_i static cache covers the whole z-range already.
    lnm_grid = np.linspace(config['lnm_low'], config['lnm_high'],
                            config['n_lnm'])
    z_grid   = np.linspace(config['zt_low'], config['zt_high'],
                            config['n_z_shared'])
    t_phod = time.perf_counter()
    lam_k, W_k, P_Mz, degenerate = _compute_lam_nodes_and_P_HOD(
        lnm_grid, z_grid, mor, config['gl_t'], config['gl_w'])
    dt_phod_ms = 1000.0 * (time.perf_counter() - t_phod)

    # K_i(lam_k, z) — closed form via CDF differencing at the unique bin
    # edges. plob_ltr_params are frozen Y3 splines, so mu/sigma/tau/fprj
    # evaluation on (lam_k, z) is pure numpy broadcasting; differencing
    # 5 edge CDFs gives the 4 unique K_i tables in one shot.
    #
    # Pass the 1-D z_grid (not a pre-broadcast (n_lnm, n_z, n_q) tensor)
    # so _plob_params can evaluate the 8 splines on 64 z-nodes only and
    # broadcast internally.  This single change saves ~130 ms/sample on
    # the (256, 64, 32) grid.
    n_lnm = lnm_grid.size
    n_z   = z_grid.size

    mu_p, sig_p, tau_p, fprj_p = _plob_params(lam_k, z_grid, plob)
    # mu_p, sig_p, tau_p, fprj_p are already broadcast to the full
    # (n_lnm, n_z, n_q) tensor via the mu = a_mu + b_mu * ltr line inside
    # _plob_params (with ltr = lam_k).
    edges = _unique_edges(config['lam_min'], config['lam_max'])   # e.g. [20,30,45,60,200]
    cdfs_at_edge = _cdf_lob_stacked(edges, mu_p, sig_p, tau_p, fprj_p)

    S_pack = np.empty((n_bins, n_z, n_lnm), dtype=np.float64)
    for k in range(n_bins):
        lo_i = int(np.searchsorted(edges, config['lam_min'][k]))
        hi_i = int(np.searchsorted(edges, config['lam_max'][k]))
        K_i = cdfs_at_edge[hi_i] - cdfs_at_edge[lo_i]
        S_i = np.sum(W_k * K_i * P_Mz, axis=-1)
        S_i = np.where(degenerate, 0.0, S_i)
        K_j_vec = _K_j(z_grid, float(config['zob_min'][k]),
                       float(config['zob_max'][k]),
                       float(config['sigma_z'][k]))
        S_pack[k] = (S_i * K_j_vec[None, :]).T

    block['sel_function', 'lnM'] = lnm_grid
    block['sel_function', 'z']   = z_grid
    block['sel_function', 'S_stack'] = np.ascontiguousarray(S_pack)

    dt_ms = 1000.0 * (time.perf_counter() - t_start)
    print(f"[sel_function] {n_bins} S_ij tables ({config['n_lnm']}x"
          f"{config['n_z_shared']}) — P_HOD {dt_phod_ms:.0f} ms, "
          f"total {dt_ms:.0f} ms", flush=True)
    return 0


def cleanup(config):
    return 0
