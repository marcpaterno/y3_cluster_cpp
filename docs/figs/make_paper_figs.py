"""Produce the three publication figures for projection_lensing_paper.tex.

Outputs (PNG, 200 dpi, bbox=tight) in this directory:
  - fig_selection_function.png   (Sec. 2.1 -- richness selection function S)
  - fig_bsel_theta.png           (Sec. 2.3 -- <b_sel,ij>(theta) probe)
  - fig_sigma_split.png          (Sec. 2.4 -- 1h+2h projected Sigma)

Notation in legends/labels matches the paper exactly:
  $\\langle b_{\\rm sel,ij}\\rangle$,
  $\\langle b_{\\rm sel,ij}^{\\rm LOS}\\rangle$,
  $\\langle b_{\\rm sel,ij}^{\\rm LSS}\\rangle$,
  $\\mathcal{S}$.

Data sources:
  - Fig 1: live richness_selection.selection_function.S_threshold.
  - Fig 2: validation TSVs at validation-data/costanzi2026/ (bsel_theta.tsv +
           bsel_plateaus.tsv) at the wall probe (lob, zob) = (25, 0.275).
  - Fig 3: live SigmaPrj for the 2h piece and NFWMiscentered.sigma_grid with
           R_mis=0 (centred) for the 1h piece, mass-weighted by P(M | lob, zob)
           inferred from the bias_precompute (lob, zob) = (20, 0.5).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# ---- locate richness_selection ----------------------------------------------
PY_SRC = "/pscratch/sd/j/jesteves/github/RichnessSelection/src"
if PY_SRC not in sys.path:
    sys.path.insert(0, PY_SRC)

from richness_selection import (
    Cosmology, PkGrid, HMF, Bias, MOR, NFWMiscentered, SelBias, SigmaPrj,
)
from richness_selection.sigma_m import SigmaM
from richness_selection.xi_nl import XiNL
from richness_selection.selection_function import S_threshold

REPO_ROOT = Path("/pscratch/sd/j/jesteves/y3_cluster_cpp")
TSV_DIR = REPO_ROOT / "validation-data" / "costanzi2026"
OUT_DIR = REPO_ROOT / "docs" / "figs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---- shared style ----------------------------------------------------------
plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "font.size": 11,
    "axes.labelsize": 11, "legend.fontsize": 9,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "xtick.top": True, "ytick.right": True,
    "legend.frameon": False, "savefig.dpi": 200, "savefig.bbox": "tight",
})
CRIMSON = "#8C1D40"


# ============================================================================
# Figure 1 -- richness selection function S(M | lambda^ob > lam_min, z=0.4)
# ============================================================================
def fig_selection_function():
    z = 0.4
    M = np.logspace(13.0, 15.3, 80)
    lam_min_list = [5, 10, 20, 30, 50]
    palette = cm.viridis(np.linspace(0.15, 0.85, len(lam_min_list)))

    mor = MOR()
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    for lam_min, col in zip(lam_min_list, palette):
        S = np.array([S_threshold(float(m), z, lam_min, mor) for m in M])
        ax.plot(M, S, color=col, lw=1.6,
                label=fr"$\lambda_{{\min}}={lam_min}$")
    ax.set_xscale("log")
    ax.set_xlim(M[0], M[-1])
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel(r"$M\;[h^{-1}M_\odot]$")
    ax.set_ylabel(r"$\mathcal{S}(M\,|\,\lambda^{\rm ob}>\lambda_{\min},\,z=0.4)$")
    ax.set_title(r"redMaPPer richness selection function (above-threshold, $z=0.4$)")
    ax.legend(loc="lower right")
    out = OUT_DIR / "fig_selection_function.png"
    fig.savefig(out)
    plt.close(fig)
    return out


# ============================================================================
# Figure 2 -- <b_sel,ij>(theta) at (lob, zob)=(25, 0.275) from validation TSVs
# ============================================================================
def fig_bsel_theta():
    """Curve = C++ ytest from bsel_theta.tsv at (lob, zob)=(25, 0.275).
    LOS / LSS asymptotes = ytest of B_small / B_large from bsel_plateaus.tsv
    at the matching (lambda_bin=0, zob_bin=0).
    """
    bsel = np.genfromtxt(TSV_DIR / "bsel_theta.tsv",
                         delimiter="\t", names=True, dtype=None, encoding="utf-8")
    plat = np.genfromtxt(TSV_DIR / "bsel_plateaus.tsv",
                         delimiter="\t", names=True, dtype=None, encoding="utf-8")

    mask = (bsel["label"] == "bsel_marg") & (bsel["lob"] == 25) & (bsel["zob"] == 0.275)
    theta_over_lam = bsel["theta_over_lam"][mask]
    y_cpp = bsel["ytest"][mask]
    order = np.argsort(theta_over_lam)
    theta_over_lam = theta_over_lam[order]
    y_cpp = y_cpp[order]

    pmask = (plat["lob"] == 25) & (plat["zob"] == 0.275)
    b_los = float(plat["ytest"][pmask & (plat["label"] == "B_small")][0])
    b_lss = float(plat["ytest"][pmask & (plat["label"] == "B_large")][0])

    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.plot(theta_over_lam, y_cpp, color=CRIMSON, lw=2.0,
            label=r"$\langle b_{\rm sel,ij}\rangle(\theta)$")
    ax.axhline(b_los, color="C0", ls="--", lw=1.4,
               label=r"$\langle b_{\rm sel,ij}^{\rm LOS}\rangle$")
    ax.axhline(b_lss, color="C2", ls=":", lw=1.4,
               label=r"$\langle b_{\rm sel,ij}^{\rm LSS}\rangle$")
    ax.axvline(0.5, color="grey", ls=":", lw=0.8)

    ax.set_xscale("log")
    ax.set_xlim(theta_over_lam.min(), theta_over_lam.max())
    ax.set_xlabel(r"$\theta/\theta_\lambda$")
    ax.set_ylabel(r"$\langle b_{\rm sel,ij}\rangle(\theta)$")
    ax.set_title(r"$\langle b_{\rm sel,ij}\rangle(\theta)$ at "
                 r"$(\lambda^{\rm ob},z^{\rm ob})=(25, 0.275)$")
    ax.legend(loc="best")
    out = OUT_DIR / "fig_bsel_theta.png"
    fig.savefig(out)
    plt.close(fig)
    return out, dict(b_los=b_los, b_lss=b_lss,
                     n_theta=len(theta_over_lam))


# ============================================================================
# Figure 3 -- 1h + 2h projected Sigma at (lob, zob) = (20, 0.5)
# ============================================================================
def fig_sigma_split():
    lob, zob = 20.0, 0.5
    R = np.logspace(-1, 1.3, 28)  # h^-1 Mpc

    cosmo = Cosmology()
    pk = PkGrid(cosmo)
    sigma_m = SigmaM(pk)
    hmf = HMF(sigma_m)
    bias = Bias(sigma_m)
    mor = MOR()
    xi_nl = XiNL(cosmo)
    sel = SelBias(cosmo, pk, hmf, bias, mor, xi_nl=xi_nl)
    nfw = NFWMiscentered(cosmo)
    sprj = SigmaPrj(cosmo, sel, nfw)

    # 2-halo (uses cl + rnd internally; the paper's <Sigma_prj> = total)
    sigma_prj = np.asarray(sprj(R, lob, zob))

    # 1-halo: centred NFW (R_mis = 0) at the bin's mass-weighted effective mass.
    # Effective mass weight: w(M) = n(M, zob) M P(lob | M, zob), normalised on ln M.
    log10_Mmin = np.log10(sel.min_mass4integral)
    m_grid = 10.0 ** np.linspace(log10_Mmin, sel.ln_M_max_log10, 200)
    n_m = hmf(m_grid, zob)
    p_lob_M = mor.pdf(np.array([lob])[:, None], m_grid[None, :], zob).ravel()
    w = n_m * m_grid * p_lob_M
    M_eff = float(np.trapz(w * m_grid, np.log(m_grid))
                  / np.trapz(w, np.log(m_grid)))

    # NFWMiscentered.sigma_grid returns shape (N_Rmis, N_R); R_mis -> 0 => centred.
    # Use a tiny positive R_mis (the table clips log to lnxmis_lo anyway).
    sigma_1h = np.asarray(nfw.sigma_grid(R, np.array([1e-6]), M_eff, zob)).ravel()

    sigma_total = sigma_1h + sigma_prj

    # Exclusion radius
    R_lambda = (lob / 100.0) ** 0.2  # cMpc/h
    R_excl = R_lambda * (1.0 + zob)

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.loglog(R, sigma_1h, color="C0", lw=2.0, ls="-",
              label=r"$\Sigma^{\rm 1h}(R)$")
    ax.loglog(R, sigma_prj, color="C2", lw=2.0, ls="--",
              label=r"$\Sigma^{\rm prj}(R)$")
    ax.loglog(R, sigma_total, color=CRIMSON, lw=2.2, ls="-",
              label=r"$\Sigma^{\rm 1h}+\Sigma^{\rm prj}$")
    ax.axvline(R_excl, color="grey", ls=":", lw=0.8)

    ax.set_xlabel(r"$R\;[h^{-1}\,\mathrm{Mpc}]$")
    ax.set_ylabel(r"$\Sigma(R)\;[h\,M_\odot/\mathrm{Mpc}^2]$")
    ax.set_title(r"Stacked surface density at "
                 r"$(\lambda^{\rm ob},z^{\rm ob})=(20, 0.5)$")
    ax.legend(loc="best")
    out = OUT_DIR / "fig_sigma_split.png"
    fig.savefig(out)
    plt.close(fig)
    return out, dict(M_eff=M_eff, R_excl=R_excl)


def main():
    p1 = fig_selection_function()
    print(f"  wrote {p1}")
    p2, info2 = fig_bsel_theta()
    print(f"  wrote {p2}  (B_LOS={info2['b_los']:.3f}, B_LSS={info2['b_lss']:.3f},"
          f" N_theta={info2['n_theta']})")
    p3, info3 = fig_sigma_split()
    print(f"  wrote {p3}  (M_eff={info3['M_eff']:.3e} h^-1 Msun,"
          f" R_excl={info3['R_excl']:.3f} h^-1 cMpc)")


if __name__ == "__main__":
    main()
