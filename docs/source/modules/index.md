# Module reference

This part documents the modules of the **reference Costanzi-2026 path**
one by one, in pipeline order. Historical, experimental, and retired
modules are confined to the {doc}`status appendix <historical>` so that no
directory in `src/modules/` is unexplained — without letting the
accumulation of past analysis choices shape the documentation.

## The template + macro module pattern

All C++ modules share one structure, explained here once and referenced
from each module section:

1. **Models** (`src/models/*.hh`, `.cuh` for device) — header-only
   templated structs, each evaluating one term of an integrand (a MOR
   probability, an HMF amplitude, an NFW $\Sigma$, a richness kernel, …).
   Each is constructed from a `cosmosis::DataBlock` and evaluated via
   `operator()`.
2. **Integrand class template** — composes the term evaluations and the
   integrator (fixed-GL evaluator, Cuhre, or PAGANI) into an object
   CosmoSIS can call.
3. **Module `.cc`** — instantiates the template and exports
   `setup`/`execute`/`cleanup` through a single
   `DEFINE_COSMOSIS_*_MODULE` macro (`src/utils/module_macros.hh`).
4. **CMake registration** — one `add_library(Name MODULE …)` per module;
   `src/modules/CMakeLists.txt` is the registry of what is built.
   (Python modules skip all of this: the ini points at the `.py` file
   directly.)

Every module section below answers, in this order: **scientific
definition**, **inputs** (DataBlock + ini knobs), **outputs**,
**numerical recipe**, **implementation**, **validation**.

## cp_camb — CosmoPower linear-P(k) emulator

### Scientific definition

`cp_camb` replaces the CAMB Boltzmann call with a CosmoPower
neural-network emulator of the linear matter power spectrum, in the
sigma8-amplitude parametrisation (the amplitude input is $\sigma_8$, not
$\ln 10^{10}A_s$). The emulators are trained at $z=0$ only; redshift
evolution is reconstructed with the linear growth factor:

$$
P(k, z) = \left[\frac{D(z)}{D(0)}\right]^2 P_{\mathrm{emu}}(k, z{=}0),
\qquad
P_{\mathrm{emu}} = 10^{\,\mathtt{NN}(h_0,\,\Omega_m,\,\Omega_b,\,n_s,\,\sigma_8,\,m_\nu)} .
$$

It also publishes background distances via astropy `FlatLambdaCDM`, so no
downstream module needs CAMB at all.

### Inputs

DataBlock reads: `cosmological_parameters/{h0, omega_m, omega_b, n_s,
sigma8, mnu}`; `growth_parameters/{z, d_z}` (when `apply_growth` is
active — requires the CosmoSIS `structure/growth_factor` module chained
*before* `cp_camb`).

| knob | production value | meaning |
|---|---|---|
| `emulator_repo` | `.../camb-emulator/camb-for-cp` | path making `cp_numpy` importable |
| `linear_pk_path` | `camb_linear_s8_v3c_emulator.npz` | total-matter linear emulator (required) |
| `linear_nonu_pk_path` | `camb_linear_nonu_s8_v3c_emulator.npz` | CDM+baryon (nonu) emulator |
| `nonlinear_pk_path` | (unset) | optional non-linear emulator |
| `nonu_fallback` | F (default) | copy linear output into `cdm_baryon_power_lin` if no nonu emulator |
| `zmin, zmax, nz` | 0.0, 4.0, 50 | output z-grid (z up to 4 so `average_sigma_crit_inv` covers the source p(z) tail) |
| `write_distances` | T (default) | publish `distances/*` via astropy |
| `apply_growth` | auto-T (z=0-only emulator) | scale $P(k,0)$ by $D(z)^2/D(0)^2$ |

### Outputs

- `matter_power_lin/{z, k_h, p_k}` — grid, `p_k` shape `(nz, n_k)`.
- `cdm_baryon_power_lin/{z, k_h, p_k}` — same grid; read by `MfTinker`
  with `matter_power_lin_version = 2`.
- `matter_power_nl/{z, k_h, p_k}` — only if a non-linear emulator is
  loaded.
- `distances/{z, a, d_a, d_m, d_l, d_c, h, mu, nz}` — Mpc (CAMB
  convention, no $h$); `d_c = d_m` (flat), `h = H(z)/c` in 1/Mpc.
- On a rejected draw: `cosmological_parameters/cp_camb_invalid_reason`
  (string) and module status 1 (→ $\log L = -\infty$).

### Numerical recipe

1. Read the 6-parameter cosmology vector; run the bound-box validator:
   reject if $\Omega_b \ge \Omega_m$ or any parameter is outside the
   emulator's trained prior (`parameters_min/max`), *before* any
   GSL-backed downstream code sees the sample.
2. Write astropy distances on the z-grid.
3. Read $D(z)$ from `growth_parameters`, renormalise by $D(0)$ (CosmoSIS
   growth is matter-dominated-normalised, $D(0)\approx 0.76$).
4. For each loaded emulator: one NN forward pass at $z=0$, then broadcast
   $P(k,z)=D(z)^2 P(k,0)$; `put_grid` the result. All emulators are
   checked at setup to share one k-grid.

### Implementation

File: `src/modules/cp_camb/cp_camb.py`. Key functions: `setup` (loads
`cp_numpy.CosmoPowerNumpyNN` artifacts, builds the bound box),
`_validate_cosmo_box`, `_read_growth`, `_evaluate`, `_write_distances`,
`execute`.

### Validation

The bound-box gate keeps the chain inside the region the emulator was
validated on and prevents GSL `qag` aborts in `cluster_toolkit` at
extrapolated cosmologies. The k-mode consistency check at setup guarantees
all published grids agree. Emulator-vs-full-CAMB residuals on the
downstream number counts are tracked in `docs/emulator_validation.tex`
(the `NumCountsFullScalarIntegrand` reference block in the ini mirrors
that report's config).

*Source: `src/modules/cp_camb/cp_camb.py`, the reference ini,
`docs/pipeline_modules.tex` §cp_camb.*

## halo_model — Tinker bias, ξ_NL, NFW lensing tables

### Scientific definition

Publishes the halo-model ingredients that half the downstream pipeline
consumes: the Tinker-2010 halo bias $b(M,z)$, the non-linear matter
correlation function $\xi_{\mathrm{NL}}(r,z)$, and (optionally) analytic
NFW one-halo and Hankel-transform two-halo lensing profiles. The bias is
computed from the $z=0$ peak height, with redshift evolution through the
growth factor:

$$
b(M, z) = b_{\mathrm{Tinker}}\!\left(\frac{\nu(M)}{D(z)/D(0)}\right),
\qquad
\nu(M) = \frac{\delta_c}{\sigma(M, z{=}0)} .
$$

### Inputs

DataBlock reads: `cosmological_parameters/{omega_m, omega_b, h0}`;
`matter_power_lin/{k_h, p_k, z}`; `matter_power_nl/*` (falls back to
linear if absent); `growth_parameters/{z, d_z}` (falls back to $D(z)=1$
if absent — wiring-test only).

| knob | production value | meaning |
|---|---|---|
| `R_perp_min, R_perp_max, R_perp_bins` | 0.05, 10.0, 128 | projected-radius grid for lensing tables (cMpc/h) |
| `Radii_min, Radii_max, Radii_bins` | 1.0, 35.0, 128 | 2h $W_p$ radius grid (cMpc/h) |
| `M_min, M_max, M_bins` | 1e12, 1e16, 100 | halo-mass grid (Msun/h) |
| `compute_lensing_1h` | T | publish NFW `Sigma_nfw/dSigma_nfw/concentration` (read by the one-halo shear module) |
| `compute_lensing_2h` | F | publish `Sigma_hh/dSigma_hh/Wp_hh` (fiducial 1h+2h composition only; skipping saves ~200–300 ms/sample) |
| `compute_lensing` | (unset) | legacy single knob, maps to both toggles |

### Outputs

- `haloModel/{m_h (100,), lnM (100,), z (50,), rhoc (50,), bias (50, 100)}`
  — `bias` is the $(z, M)$ table read via `Interp2D` by `b_sel_marg` and
  `shear_prj`.
- `xi_nl/{r (128,), z (50,), xi_nl (50, 128)}` — $\xi_{\mathrm{NL}}$ on
  $r \in [10^{-3}, 10^{3}]$ cMpc/h, always written.
- If 1h: `haloModel/{r_sigma, scale_shift, hubble_shift, k, Sigma_nfw,
  dSigma_nfw, concentration}`.
- If 2h: `haloModel/{Rp, Wp_hh, Sigma_hh, dSigma_hh}`.

### Numerical recipe

1. Disable GSL's abort-on-error handler at import (via ctypes) so extreme
   cosmologies raise Python exceptions instead of killing the worker.
2. $\nu(M)$ once at $z=0$ from `cluster_toolkit.peak_height.nu_at_M` on
   the linear $P(k)$; per z-slice, evaluate
   `bias_at_nu(nu / (D(z)/D(0)))` ($\Delta = 200$).
3. Per z-slice, $\xi_{\mathrm{NL}}(r,z)$ from `ct.xi.xi_mm_at_r` on the
   (non-linear, or fallback linear) $P(k)$.
4. Optional lensing branches: `lensingModel.first_halo_term` (Child-18
   concentration NFW $\Sigma/\Delta\Sigma$) and `second_halo_term`
   (Hankel $P \to \xi \to \Sigma \to \Delta\Sigma$).

### Implementation

File: `y3_buzzard/halo_model_cosmosis.py`. Uses
`haloModel.{biasModel, lensingModel, scaleShiftCosmo}` and
`cluster_toolkit`; entry points `setup`/`execute`.

### Validation

The $D(0)$ renormalisation is load-bearing: pre-fix, using un-normalised
$D(z)$ inflated $\nu$ by $1/D(0)\simeq 1.32$ and produced halo biases up
to 2× too large (documented in `docs/pipeline_modules.tex`). Bias and
$\xi_{\mathrm{NL}}$ tables are cross-checked against the Python
`richness_selection` reference through the `compare_*_py_vs_cpp.py`
harnesses in `$PSCRATCH/github/RichnessSelection/validations/`.

*Source: `y3_buzzard/halo_model_cosmosis.py`, the reference ini,
`docs/pipeline_modules.tex` §halo_model.*

## average_sigma_crit_inv — source-averaged $\Sigma_{\mathrm{crit}}^{-1}$

### Scientific definition

Computes the source-population-weighted mean inverse critical surface
density as a function of lens redshift, the geometric factor that
converts $\Delta\Sigma$ into tangential shear:

$$
\langle \Sigma_{\mathrm{crit}}^{-1} \rangle(z_l)
= h_0 \int dz_s\; p(z_s + \delta_z)\,
  \frac{4\pi G}{c^2}\,
  \frac{D_A(z_l)\, \big[ D_A(z_s) - \frac{1+z_l}{1+z_s} D_A(z_l) \big]}{D_A(z_s)} ,
$$

clipped at zero for sources in front of the lens. With `unity = T` the
module instead publishes $\langle \Sigma_{\mathrm{crit}}^{-1} \rangle
\equiv 1$, so the downstream shear modules
($\gamma_t = \Delta\Sigma \cdot \Sigma_{\mathrm{crit}}^{-1}$) emit
$\Delta\Sigma$ directly — used when the mock data vector is a Buzzard
jackknife $\Delta\Sigma$ measurement.

### Inputs

DataBlock reads: `cosmological_parameters/h0`; `photoz/delta_z` (source
photo-z shift nuisance); `distances/{z, d_a}` (Mpc, from `cp_camb`). The
source $n(z)$ is read at setup from
`${Y3_CLUSTER_CPP_DIR}/data/test_cluster_Y1.fits` (HDU 6, columns
`z_mid`, `bin1`).

| knob | production value | meaning |
|---|---|---|
| `z_min, z_max` | 0.05, 0.80 | lens-redshift grid range |
| `z_bins` | 50 | number of lens-z points |
| `unity` | T (production mock closure) | publish `sci_average = 1` (ΔΣ-observable mode) |

### Outputs

- `average_sigma_crit_inv/zlense` — `(z_bins,)` lens-redshift grid.
- `average_sigma_crit_inv/sci_average` — `(z_bins,)`
  $\langle \Sigma_{\mathrm{crit}}^{-1} \rangle(z_l)$ (or ones if
  `unity = T`).

### Numerical recipe

For each $z_l$ on the linear lens grid: interpolate $D_A$ (1-D `interp1d`
on `distances`), shift the source p(z) by `delta_z`, form the integrand
$\max(0,\, p(z_s{+}\delta_z)\, 4\pi G c^{-2} D_A(z_l) D_A(z_l, z_s)/D_A(z_s))$,
and integrate with `np.trapz` over the source grid; multiply by $h_0$.
Constants are in Mpc/Msun/s units
($G = 4.517\times 10^{-48}\,\mathrm{Mpc^3\,M_\odot^{-1}\,s^{-2}}$).

### Implementation

File: `src/modules/average_sigma_crit_inv/average_sigma_crit_inv.py`.
Plain `setup`/`execute`; the `unity` branch short-circuits `execute`
before any distance work.

### Validation

`execute` prints $\sum p(z_s)$ as a normalisation sanity check on the
FITS n(z). The `unity` mode is the closure switch for the $\Delta\Sigma$
mock-DV pipeline: with it on, the shear observable is bit-comparable to
the jackknife $\Delta\Sigma$ data with no lensing-geometry systematics in
play.

*Source: `src/modules/average_sigma_crit_inv/average_sigma_crit_inv.py`,
the reference ini.*

## sel_function — $S_{ij}(\ln M, z)$ tabulation

### Scientific definition

Pre-tabulates, once per sample, the richness–photo-z selection function
for all 12 wall bins on one shared $(\ln M, z)$ grid; the C++
`NumCountsSel`/one-halo shear integrands slice their bin's plane and
serve it via `Interp2D`:

$$
\begin{aligned}
S_{ij}(\ln M, z) &= \Big[ \textstyle\sum_k W_k\, K_i(\lambda_k, z)\,
  P_{\mathrm{HOD}}(\lambda_k \mid M, z) \Big] \cdot K_j(z), \\
K_i(\lambda^{\mathrm{tr}}, z) &= F_{\mathrm{EMG}}(\lambda_{\max}^i \mid \lambda^{\mathrm{tr}}, z)
  - F_{\mathrm{EMG}}(\lambda_{\min}^i \mid \lambda^{\mathrm{tr}}, z), \\
K_j(z) &= \Phi\!\big((z_{\max}^j - z)/\sigma_z\big) - \Phi\!\big((z_{\min}^j - z)/\sigma_z\big),
\end{aligned}
$$

with $P_{\mathrm{HOD}}$ the continuous shifted-Poisson form
$\exp[-\nu + (\lambda^{\mathrm{tr}} + \delta - 1)\ln\nu -
\ln\Gamma(\lambda^{\mathrm{tr}} + \delta)]$,
$\nu = \mu_{\mathrm{sat}} + \delta$,
$\delta = (\sigma_\lambda \mu_{\mathrm{sat}})^2$, matching
`src/models/mor_hod_t.hh` line-for-line.

### Inputs

DataBlock reads: `cluster_mor/{log10_Mmin, log10_M1 | log10_ratio,
alpha, epsilon, sigma_lambda, z_pivot?}`; EMG coefficient splines from
`plob_ltr_params/*` when present, with fallback to the in-code
`y3_buzzard.prj_params.PrjParams` table.

| knob | production value | meaning |
|---|---|---|
| `lam_min, lam_max` | {20,30,45,60} → {30,45,60,200} ×3 z-bins | per-bin richness edges (12 bins) |
| `zob_min, zob_max` | {0.20,0.35,0.50} → {0.35,0.50,0.65} | per-bin photo-z edges |
| `sigma_z` | 0.03 (all bins) | photo-z scatter in $K_j$ |
| `zt_low, zt_high` | 0.05, 0.80 | shared z-grid envelope |
| `lnm_low, lnm_high` | 29.9336, 36.8414 | shared lnM-grid envelope |
| `n_lnm` | 192 | lnM nodes (whole-pipeline optimum; 64 is pathological — GL resonance, 4.5% drift) |
| `n_z_shared` | 64 | shared z nodes |
| `N_q` | 32 | GL nodes for the $\lambda^{\mathrm{tr}}$ quadrature |
| `L_lam`, `L_z` | 6.0, 6.0 | bracket half-widths in units of $\sigma$ |

### Outputs

- `sel_function/lnM` — `(192,)`; `sel_function/z` — `(64,)`.
- `sel_function/S_stack` — `(12, 64, 192)` packed tensor,
  `(bin, z, lnM)`, C-contiguous.

### Numerical recipe

1. Per $(\ln M, z)$ cell, place $N_q$ Gauss–Legendre nodes in
   $\lambda^{\mathrm{tr}}$ on
   $[\mu_{\mathrm{eff}} - L\sigma_{\mathrm{eff}},\,
   \mu_{\mathrm{eff}} + L\sigma_{\mathrm{eff}}]$,
   $\sigma_{\mathrm{eff}} = \sqrt{\mu_{\mathrm{sat}} +
   (\sigma_\lambda\mu_{\mathrm{sat}})^2}$; evaluate $P_{\mathrm{HOD}}$ on
   the full `(192, 64, 32)` tensor in a single `gammaln` call
   (narrow-Gaussian fallback where $\mu_{\mathrm{sat}} \le 10^{-8}$).
2. Evaluate the 8 EMG coefficient splines on the 1-D z-grid only (fast
   path, saves ~130 ms/sample), broadcast to
   $(\mu, \sigma, \tau, f_{\mathrm{prj}})$ at each node.
3. Compute the EMG CDF via `erfcx` at the 5 unique bin edges
   {20, 30, 45, 60, 200} and difference to get all 4 $K_i$ tables at
   once.
4. Per bin: contract $S_i = \sum_k W_k K_i P_{\mathrm{HOD}}$, multiply by
   $K_j(z)$, pack into `S_stack`.

### Implementation

File: `src/modules/sel_function/sel_function.py`. Key functions:
`_compute_lam_nodes_and_P_HOD`, `_p_hod_scalar`, `_plob_params`,
`_cdf_lob` / `_cdf_lob_stacked`, `_K_edges_of_bins`, `_K_j`, `execute`.

### Validation

`_f_emg`/`_cdf_lob` match `src/models/richness_kernel_t.hh::F_EMG`;
`_p_hod_scalar` matches `mor_hod_t.hh`; `_K_i_bin` is kept as a per-bin
parity/debug path bit-identical to the stacked production path. The
`n_lnm` sweep (2026-05-07) shows <0.05% change on NumCounts and the
one-halo shear at 192 vs 256; `execute` prints per-sample timing.

*Source: `src/modules/sel_function/sel_function.py`, the reference ini,
`docs/pipeline_modules.tex` §sel_function.*

## NumCountsSel — cluster number counts

### Scientific definition

`NumCountsSel` computes the expected cluster number count $N_i[1]$ in
each of the 12 DES Y3 $(\lambda^{\rm ob}, z^{\rm ob})$ bins, using the
pre-tabulated richness-selection tensor
$S_{ij}(\ln M, z) = S_i(\ln M, z)\,K_j(z)$ published by `sel_function`.
It is one instantiation of the generic $N_i[f]$ operator (with weight
$f = 1$):

$$
N_i[f] = \int d\ln M \int dz\;
  \Omega(z)\,\frac{dV}{d\Omega\,dz}\,n(M, z)\,
  S_{ij}(\ln M, z)\, f(\ln M, z).
$$

The full derivation lives in {doc}`../science/index`; here only the
operator contract matters.

### Inputs

DataBlock sections read (in `set_sample`):

- `sel_function/{lnM, z, S_stack}` — per-bin $S_{ij}$ tables (bin count =
  first extent of `S_stack`);
- `mass_function/{m_h, z, dndlnmh}` — via `HMF_t`, which rescales the
  mass axis by $(\Omega_m - \Omega_\nu)$ internally (so `ln_mass_shift`
  must stay 0);
- `cluster_abundance/{hmf_s, hmf_q}` — HMF nuisance amplitudes;
- `distances/{z, d_a}` — via `DV_DO_DZ_t` for $dV/d\Omega\,dz$;
- `cosmological_parameters/{omega_M, omega_nu}`.

| knob | production value | meaning |
|---|---|---|
| `bin_index` | `0 1 ... 11` | wall of bins (richness-fast, $z$-block-slow) |
| `zt_low`, `zt_high` | 0.05, 0.80 | true-$z$ integration limits |
| `lnm_low`, `lnm_high` | 29.9336, 36.7300 | $\ln M$ integration limits |
| `n_lnm` | 96 (default) | GL nodes in $\ln M$ |
| `n_z` | 64 (default) | GL nodes in $z$ ($S_{ij}$ has compact $\sim 0.15$-wide $z$ support) |
| `algorithm`, `eps_rel`, `eps_abs`, `max_eval`, `use_cartesian_product` | (present in ini) | legacy Cuhre knobs, **ignored** by the fixed-GL evaluator |

### Outputs

- `numcountssel/vals` — length-12 array, $N_i[1]$ per bin. The output
  section is hard-coded (deliberately *not* an ini knob: CosmoSIS
  `[DEFAULT]` blocks would propagate `output_section` into every module
  and silently redirect the write).

### Numerical recipe

Fixed Gauss–Legendre replaces the retired per-bin adaptive Cuhre
integral. Once per sample, `SelGLCore::build_weights` contracts the $z$
axis:

$$
W_{ij}(\ln M_k) = \sum_q w_q\,
  \frac{dV}{d\Omega\,dz}(z_q)\,\Omega(z_q)\,
  n(M_k, z_q)\, S_{ij}(\ln M_k, z_q),
$$

on `n_z` GL nodes; each bin's count is then the 1-D GL sum
$N_i = \sum_k w_k\, W_{ij}(\ln M_k)$. Cost is deterministic (0.021 s vs
Cuhre mean 0.107 s / max 0.98 s over $\sim 10^6$ MCMC realisations).
Avoid `sel_function/n_lnm = 64`, which hits a GL resonance ($-4.55\%$ on
$N_c$); the whole-pipeline optimum is 192.

### Implementation

- `src/models/n_operator_sel_gl_t.hh` — `y3_cluster::NumCountsSelGL`
  wrapping `nosel_gl_detail::SelGLCore`;
- `src/modules/num_counts_sel/NumCounts.cc` —
  `DEFINE_COSMOSIS_SCALAR_EVALUATOR_MODULE(NumCountsSel)`;
- sibling instantiations of the older Cuhre-driven template
  `NOperatorSelScalar<F>` (`src/models/n_operator_sel_t.hh`):
  **MassWeightedSel** (`MassWeighted.cc`, $f = M$, gives
  $\langle M\rangle_i = N_i[M]/N_i$) and **BiasWeightedSel**
  (`BiasWeighted.cc`, $f = b(M,z)$ read from `haloModel/{lnM, z, bias}`,
  gives $\langle b\rangle_i$), both via
  `DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE`;
- brute-force reference:
  `src/modules/num_counts_full/NumCountsFullScalarIntegrand.cc`
  integrates the full triple $(\lambda^{\rm tr}, z, \ln M)$ integral with
  `RmSelFunction_t` $B_i$ kernels and Cuhre (no miscentering, no
  photo-$z$), one cell per bin of the 12-cell wall.

### Validation

`NumCountsFullScalarIntegrand` is the brute-force cross-check for
anything landing on the same datablock layout (its config mirrors the
emulator-validation report so cp_camb-vs-CAMB residuals match
`docs/emulator_validation.tex` §2). The `n_lnm` sweep bounds the GL-grid
error at $< 0.05\%$ (192 vs 256 nodes). Python↔C++ number-count
harnesses live in `$PSCRATCH/github/RichnessSelection/validations/`.

*Source: `src/models/n_operator_sel_gl_t.hh`,
`src/models/n_operator_sel_t.hh`, `src/modules/num_counts_sel/`,
`src/modules/num_counts_full/`, the reference ini,
`docs/pipeline_modules.tex`.*

## Shear1hMisSel — one-halo shear with miscentering

### Scientific definition

`Shear1hMisSel` is the production 1-halo shear branch: the $N_i[f]$
operator with the centred + miscentred tangential-shear weight. It is
*population-integrated*, not per-cluster — the likelihood divides by
`numcountssel/vals` to form
$\langle\gamma_t^{1h}\rangle_i(R) = N_i[\gamma_t^{1h,\rm full}](R) / N_i[1]$,
then adds $\gamma_t^{\rm prj}$ from `shear_prj`:

$$
\begin{aligned}
\gamma_t^{1h,\rm full}(R; M, z) ={}&
  \Big[(1 - f_{\rm mis})\,\Delta\Sigma_{\rm NFW}(R, M) \\
  &{}+ f_{\rm mis}\,\Delta\Sigma_{\rm mis}\big(R, M;\, \tau_{\rm mis} R_\lambda\big)\Big]\,
  \Sigma_{\rm crit}^{-1}(z),
\end{aligned}
$$

with the gamma-kernel (Kelly et al. 2023) offset distribution for
$\Delta\Sigma_{\rm mis}$ and
$R_\lambda = (\lambda/100)^{0.2}\,h^{-1}\mathrm{Mpc}$. Both pieces are
linear in $\Delta\Sigma$ (the reduced-shear denominator was retired
2026-05-11), so the 1h + prj sum in the likelihood is exact.

### Inputs

DataBlock: everything `NumCountsSel` reads, plus

- `haloModel/{r_sigma, lnM, dSigma_nfw}` — centred NFW $\Delta\Sigma$
  spline (requires `halo_model` run with `compute_lensing_1h = T`);
- `average_sigma_crit_inv/{zlense, sci_average}` —
  $\langle\Sigma_{\rm crit}^{-1}\rangle(z)$, folded into the $z$ weight;
- `miscentering/f_mis`, `miscentering/tau_mis` — mixture parameters,
  in-code defaults 0.22 / 0.17 if absent;
- `cosmological_parameters/omega_M` — sets the $\rho_{\rm mean}$
  multiplier of the miscentred table;
- disk tables `data/nfw_off_center/*gamma*` ($1000\times1000$ log-log
  grids in $(R/r_s, R_{\rm mis}/r_s)$, loaded once at module
  construction).

| knob | production value | meaning |
|---|---|---|
| `bin_index` × `r_perp` | 0–11 × 10 radii | Cartesian wall grid, bin slow / $R$ fast |
| `zt_low`, `zt_high` | 0.05, 0.80 | $z$ limits |
| `lnm_low`, `lnm_high` | 29.9336, 36.7300 | $\ln M$ limits |
| `n_lnm`, `n_z` | 96, 64 (defaults) | GL nodes |
| `lob_centers` | 25 37.5 52.5 130 (default) | drives $R_\lambda$ per richness bin |
| `method` | `exact` (default) | `exact` = full GL mass sum; `idea2` = 2nd-order moment expansion |
| `stencil_h` | 0.15 (default) | $\ln M$ stencil for the `idea2` second derivative |

### Outputs

- `shear1hmissel/vals` — length 120 ($12 \times 10\,R$),
  $N_i[\gamma_t^{1h,\rm full}](R)$. Section hard-coded.

### Numerical recipe

The profile $\Phi_i(R, \ln M)$ is $z$-free, so the $z$-marginalised
weight $W_{ij}(\ln M)$ — now including the $\Sigma_{\rm crit}^{-1}(z)$
factor — is built once per sample on fixed GL nodes; each of the 120 grid
points is then a single 1-D GL mass sum
$\sum_k w_k\, W_{ij}(\ln M_k)\,\Phi_i(R, \ln M_k)$. This replaces one
adaptive Cuhre integral per (bin, $R$) pair ($\sim 16\times$ faster,
deterministic cost; Cuhre mean 0.575 s, max 4.0 s). The richness bin for
$R_{\rm mis} = \tau_{\rm mis} R_\lambda$ is `bin_index % 4` — a
deliberate fix; the old weight silently reused bin 3's $R_\lambda$ for
bins 4–11 (changes those bins by up to $\sim 2\%$).

### Implementation

- `src/models/n_operator_sel_gl_t.hh` — `y3_cluster::Shear1hMisSelGL`
  (shares `SelGLCore` with `NumCountsSelGL`; miscentred table via
  `NFW_DSIGMA_MIS` with $c = 4$ and the gamma kernel);
- `src/modules/num_counts_sel/Shear1hMis.cc` —
  `DEFINE_COSMOSIS_SCALAR_EVALUATOR_MODULE(Shear1hMisSel)`;
- centred-only variant **Shear1hSel** (`Shear1h.cc`):
  `NOperatorSelRadial<Shear1hWeight>` from
  `src/modules/num_counts_sel/lensing_weights.hh`, Cuhre-driven via
  `DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE`, writes `shear1hsel/vals`.
  The mock data vector was generated against this branch.

### Validation

Setting `miscentering/f_mis = 0` recovers the centred `Shear1hSel`
result, closing the loop against the pre-baked mock. The radial
factorisation and its error budget are documented in
`docs/shear1h_radial_factorization.tex`; at the fiducial
$(f_{\rm mis}, \tau_{\rm mis}) = (0.22, 0.17)$ the small-$R$ signal is
suppressed $\sim 30\%$ at $R \lesssim 0.3\,h^{-1}$Mpc, matching the
Costanzi-2026 reference figure.

*Source: `src/models/n_operator_sel_gl_t.hh`,
`src/modules/num_counts_sel/{Shear1hMis,Shear1h}.cc`,
`src/modules/num_counts_sel/lensing_weights.hh`,
`docs/pipeline_modules.tex` §Miscentering selection,
`docs/shear1h_radial_factorization.tex`, the reference ini.*

## b_sel_marg — the $P[X]$ operators

### Scientific definition

`b_sel_marg` co-computes the three Costanzi-2026 $P[X]$ scalars on the
12-bin $(z^{\rm ob}, \lambda^{\rm ob})$ wall, feeding the analytic
`bsel` closure for the selection-bias plateaus
$(B_{\rm small}, B_{\rm large})$:

$$
\begin{aligned}
P_1 &= \int dz\, d\ln M\, d\lambda^{\rm tr}\, d\theta\;
       \mathcal{W}\, f_A\big(\theta, \theta_\lambda(\lambda^{\rm tr})\big), \\
I_1 &= \int dz\, d\ln M\, d\lambda^{\rm tr}\, d\theta\;
       \mathcal{W}\, b(M,z)\,\xi_{\rm NL}(|\Delta\chi|, z^{\rm ob})\,
       \sigma(\theta)\, f_A, \\
J   &= I_2 - I_1
     = \int dz\, d\ln M\, d\lambda^{\rm tr}\, d\theta\;
       \mathcal{W}\, b\,\xi_{\rm NL}\,\big(1 - \sigma(\theta)\big)\, f_A,
\end{aligned}
$$

with $\mathcal{W} = w_z\, (dV/d\Omega\,dz)\,
\mathbb{1}[\theta > \theta_{\rm excl}(z)]\, (dn/d\ln M)\,
P_{\rm HOD}(\lambda^{\rm tr}|M,z)\, \lambda^{\rm tr}\, 2\pi\sin\theta$
and the sigmoid $\sigma(\theta) = [1 + e^{-k(\theta - \theta_0)}]^{-1}$,
$k = 2.5/\theta_\lambda$, $\theta_0 = \theta_\lambda/2$.

### Inputs

DataBlock: `distances/{z, d_c}` (times $h_0$ for cMpc/$h$),
`mass_function/*` (via `HMF_t`), `haloModel/{lnM, z, bias}`,
`xi_nl/{r, z, xi_nl}`, `cluster_mor/*` (shifted-Poisson HOD parameters
via `MOR_HOD_t`: `log10_Mmin`, `log10_M1`, `alpha`, `epsilon`,
`z_pivot`, `sigma_lambda`), `cosmological_parameters/{omega_M, omega_nu,
h0}`.

| knob | production value | meaning |
|---|---|---|
| `zo_low`, `zo_high`, `lambda_bin` | 12-entry wall | grid; $z^{\rm ob} = (z_{\rm lo}+z_{\rm hi})/2$, richness centre from `lob_center(bin)` |
| `lnm_low`, `lnm_high` | 29.9336, 35.6814 | $\ln M$ GL limits |
| `n_lt` | 60 | GL nodes in $\lambda^{\rm tr}$ on $(0, \lambda^{\rm ob}_{\rm centre}]$ per bin (converged at 16) |
| `n_lnm` | 24 | GL nodes in $\ln M$ |
| `n_theta` | 10 | $\theta$ GL nodes, split at $\theta_\lambda$ |
| `n_zring` | 20 | GL nodes on the ring band $[z^{\rm ob} \pm \Delta z_{\rm excl}]$ |
| `n_zouter` | 20 (default 30) | GL nodes per fg/bg $\log|\Delta\chi|$ wing |
| `n_zt_ref`, `zt_ref_low`, `zt_ref_high` | 80, 0.05, 0.90 (defaults) | reference $z$ grid for per-sample tensor caches |

### Outputs

Three sections, each a length-12 `vals` array: `b_sel_marg_P1/vals`,
`b_sel_marg_I1/vals`, `b_sel_marg_J/vals`. $J$ is emitted directly (not
as $I_2 - I_1$) because $\sigma(\theta) \to 1$ at large $\theta$ makes
the difference catastrophically cancel — and $J$ is the denominator of
$b_{\rm zero}$ downstream in `bsel.py`.

### Numerical recipe

$\theta$-outer fixed GL. Per sample (`set_sample`): pre-evaluate the
shifted-Poisson MOR pdf on a $(4, N_{\rm lt}, N_{\ln M}, N_{z,\rm ref})$
tensor plus HMF/bias/$\chi$/$dV$/$\sigma_z$ reference tables; copy
$\xi_{\rm NL}$ into a flat table for a hand-rolled log-$r$ bilinear
lookup ($\sim 10\times$ faster than GSL). Per grid point: (i) build the
$z$ grid as ring (GL in $z$ on $[z^{\rm ob} \pm \Delta z_{\rm excl}]$) +
fg/bg wings (GL in $u = \ln|\Delta\chi_\parallel|$ with $du \to dz$
Jacobian), photo-$z$ endpoints by bisection of
$z \pm \sigma_z(z) = z^{\rm ob}$; (ii) build ONE $\theta$ grid on
$[\min_z\theta_{\rm excl}(z), 2\theta_\lambda]$ split at
$\theta_\lambda$, caching $\sin$, $\cos$, $\sigma(\theta)$,
$w_\theta 2\pi\sin\theta$ once; (iii) pre-contract the mass integral into
$(N_{\rm lt}, N_z)$ caches $W_{P_1}, W_I$; (iv) run $\theta$-outer /
$z$-inner / $\lambda^{\rm tr}$-innermost multiply-adds with per-$z$
exclusion gating $\theta > \theta_{\rm excl}(z)$,
$\cos\theta_{\rm excl} = (\chi_z^2 + \chi_o^2 - R_{\rm excl}^2)/(2\chi_z\chi_o)$.
`ln_mass_shift_ = 0.0` always (`HMF_t` rescales internally). $\sim 74$ ms
for all 12 bins.

### Implementation

- `src/models/p_operator_t.hh` — `y3_cluster::P_operator`
  (`module_label() = "b_sel_marg"`, `n_outputs = 3`);
- `src/modules/b_sel_marg_cpu/BSelMargIntegrand.cc` —
  `DEFINE_COSMOSIS_SCALAR_EVALUATOR_MODULE(BSelMargIntegrand)`;
- GPU/adaptive reference benchmarks in the same directory:
  `P1PaganiIntegrand.cc`, `I1PaganiIntegrand.cc`, `I2PaganiIntegrand.cc`
  — `P_operator_cuhre<BSelMarg*PaganiWeight>` via
  `DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE`, one adaptive integral per
  (bin, operator). Reference-only: the fixed-GL path is
  $\sim 10^3\times$ faster on this wall grid.

### Validation

Mirrors `richness_selection/sel_bias.py::_P_operator`; post-audit
agreement $(P_1, I_1, J) \leq 1\%$ against the Python reference. The
quadrature recipe (ring/wing splits, node counts, sub-0.01% convergence
vs `scipy.integrate.quad`) follows
`RichnessSelection/docs/richness_selection.tex`; the sweep harness is
`sweep_p_operator.py` in `$PSCRATCH/github/RichnessSelection/validations/`.

*Source: `src/models/p_operator_t.hh`,
`src/modules/b_sel_marg_cpu/`, `docs/pipeline_modules.tex` §b_sel_marg,
the reference ini.*

## bsel — analytic sigmoid closure for $b_{\rm sel}(\theta)$

### Scientific definition

Python finisher of the Costanzi-2026 selection-bias pipeline: it reads
the three co-computed operators $(P_1, I_1, J{=}I_2{-}I_1)$ from the C++
`BSelMargIntegrand`, closes the small- and large-scale bias asymptotes
per latent richness, marginalises over $\lambda^{\mathrm{tr}}$, and
publishes two scalars per bin from which the full scale dependence is
analytic:

$$
\begin{aligned}
\Delta_{\mathrm{RND}} &= P_1 + b_{\mathrm{eff}} I_2, \qquad
\delta_{\mathrm{prj}}(\lambda^{\mathrm{tr}}) =
  \frac{\lambda^{\mathrm{ob}} - \lambda^{\mathrm{tr}}}{\Delta_{\mathrm{RND}}} - 1, \\
b_\infty &= b_{\mathrm{eff}}\,(1 + 0.13\,\delta_{\mathrm{prj}}), \qquad
b_{\mathrm{zero}} =
  \frac{(\lambda^{\mathrm{ob}} - \lambda^{\mathrm{tr}}) - P_1 - b_\infty I_1}{J}, \\
\langle b_{\mathrm{sel}}\rangle(\theta) &= B_{\mathrm{small}} +
  (B_{\mathrm{large}} - B_{\mathrm{small}})\,\sigma(\theta), \qquad
\sigma(\theta) = \big[1 + e^{-k(\theta - \theta_0)}\big]^{-1},
\end{aligned}
$$

with $k = 2.5/\theta_\lambda$, $\theta_0 = \theta_\lambda/2$,
$\theta_\lambda = R_\lambda(\lambda^{\mathrm{ob}})(1+z^{\mathrm{ob}})/\chi(z^{\mathrm{ob}})$,
and $B_{\mathrm{small}/\mathrm{large}} =
\langle b_{\mathrm{zero}/\infty}\rangle_{\lambda^{\mathrm{tr}}}$. Since
only the sigmoid depends on $\theta$, the $\lambda^{\mathrm{tr}}$
marginalisation factorises exactly into these two scalars — no $\theta$
tabulation is needed downstream.

### Inputs

DataBlock reads: `b_sel_marg_P1/vals`, `b_sel_marg_I1/vals`,
`b_sel_marg_J/vals` (each flat `(N_zob·N_lob,)`, λ-bin index fastest);
`mass_function/{m_h, z, dndlnmh}` (mass axis rescaled by
$\Omega_m - \Omega_\nu$ to physical M/h); `halomodel/{m_h, z, bias}`;
`cluster_mor/{log10_Mmin, log10_M1 | log10_ratio, alpha, epsilon,
sigma_lambda}`; `distances/{z, d_c}` and `cosmological_parameters/h0`
(χ converted to cMpc/h); EMG coefficients from `PrjParams.default()`.

| knob | production value | meaning |
|---|---|---|
| `lob` | 25.0 37.5 52.5 130.0 | richness bin centres |
| `zob` | 0.275 0.425 0.575 | photo-z bin centres |
| `n_theta, theta_lo, theta_hi` | 32, 1e-4, 5e-3 | legacy tabulated θ grid (rad, geometric) |
| `n_ltr, ltr_lo, ltr_hi_factor` | 128, 1.0, 3.0 | GL marginalisation on $[1,\, 3\,\lambda^{\mathrm{ob}}]$ per bin |
| `ltr_hi` | 0 (off) | fixed upper edge override |
| `min_mass4integral, ln_M_max_log10, n_m_beff` | 1e13, 15.5, 100 | b_eff mass grid |
| `verbose` | F | timing + asymptote tables |

### Outputs

Section `b_sel_marginalised`: `lob (4,)`, `zob (3,)`, `theta (32,)`,
`vals (4, 3, 32)` (backward-compat tabulation), `b_eff (3, 4)`,
`b_small (3, 4)`, `b_large (3, 4)`. The `shear_prj` C++ evaluator
consumes the two-scalar form and rebuilds $b_{\mathrm{sel}}(\theta)$
analytically.

### Numerical recipe

Per $(z^{\mathrm{ob}}, \lambda^{\mathrm{ob}})$ bin: (1)
bilinear-interpolate $n(M,z)$ and $b(M,z)$ onto a 100-point log-M grid
($dndlnmh$ divided by $M$ to get $dn/dM$); (2)
$b_{\mathrm{eff}} = \int dM\, n\,
P(\lambda^{\mathrm{tr}}{=}\lambda^{\mathrm{ob}}|M,z)\, M\, b \,/\,
\int dM\, n\, P\, M$, trapezoidal in $\ln M$; (3) build 128 GL nodes in
$\lambda^{\mathrm{tr}}$, evaluate $b_{\mathrm{zero}}, b_\infty$
closed-form using $J$ directly as the denominator (C++ emits $J$ to
dodge $I_2 - I_1$ cancellation at $\sigma\to 1$; a collapsed $J$
surfaces NaN rather than silently folding $b_{\mathrm{LOS}}$ into
$b_{\mathrm{LSS}}$); (4) marginalisation weight
$w = W_{\mathrm{GL}} \cdot
P_{\mathrm{EMG}}(\lambda^{\mathrm{ob}}|\lambda^{\mathrm{tr}}, z) \cdot
\int dM\, n\, P_{\mathrm{HOD}}(\lambda^{\mathrm{tr}}|M,z)\, M$.

### Implementation

File: `y3_buzzard/bsel.py`. Key functions: `_compute_b_eff`,
`_ltr_weight`, `_p_lob_given_ltr_emg`, `_p_ltr_given_M`, `_sigmoid`,
`_theta_lob`, `execute`.

### Validation

Matches the Python reference
`richness_selection.sel_bias.SelBias.b_sel_marginalised` (with
`use_plob_ltr=True`). Two historical bugs are guarded in-code: the $h_0$
factor on $\chi$ (raw Mpc shifted the sigmoid by $h_0$ in θ, a 12%
residual at the lob=25, zob=0.275 probe, fixed 2026-05-19) and the
$dn/dM$ vs $dn/d\ln M$ convention (extra factor of $M$ halved
$B_{\mathrm{small}}$, fixed 2026-05-06). `verbose = T` prints the
$b_{\mathrm{eff}}/B_{\mathrm{small}}/B_{\mathrm{large}}$ tables per bin
for smoke runs.

*Source: `y3_buzzard/bsel.py`, the reference ini,
`docs/pipeline_modules.tex` §bsel.*

## shear_prj — two-halo projection lensing

### Scientific definition

`shear_prj` (class `ShearPrjEvaluator`, module label `shear_prj`) is the
Costanzi-2026 two-halo projection module: it evaluates, per
$(\lambda^{\rm ob}, z^{\rm ob}, R)$ wall point, the projected surface
density

$$
\begin{aligned}
\Sigma^{\rm prj}(R \,|\, \lambda^{\rm ob}, z^{\rm ob})
  = \int dz\, d\ln M\, d\theta\;&
    w_z(z, z^{\rm ob})\, \frac{dV}{d\Omega\,dz}\, n(M, z)\,
    \big[1 + b(M,z)\, b_{\rm sel}(\theta)\, \xi_{\rm NL}\big] \\
  &\times \Sigma_{\rm mis}\big(R \,|\, M, z,\, R_{\rm mis} = \theta D_A(z^{\rm ob})\big)\,
    \mathbb{1}\big[\theta > \theta_{\rm excl}(z)\big],
\end{aligned}
$$

together with $\Delta\Sigma^{\rm prj}$ and the tangential shear
$\gamma_t^{\rm prj} = \Delta\Sigma^{\rm prj}\,\langle\Sigma_{\rm crit}^{-1}\rangle$
(linear — the reduced-shear denominator was retired 2026-05-11). Unlike
the $N_i[f]$ family it uses a parabolic photo-$z$ kernel
$w_z = \max(0, 1-u^2)$ instead of $S_{ij}$, applies no $\Omega(z)$
factor, and uses the **single**-offset (delta-kernel) miscentred NFW —
here the offset *is* the $\theta$ integration variable, so a gamma
kernel would double-integrate.

### Inputs

DataBlock: `haloModel/{lnM, z, bias}`, `xi_nl/{r, z, xi_nl}`,
`distances/{z, d_c}`, `mass_function/*` (via `HMF_t`),
`b_sel_marginalised/{lob, zob, theta, vals, b_small, b_large}` (written
by `bsel.py`; the $(B_{\rm small}, B_{\rm large})$ scalars are linearly
interpolated to the exact $z^{\rm ob}$ and $b_{\rm sel}(\theta)$ is
evaluated analytically), `average_sigma_crit_inv/{zlense, sci_average}`
(optional — if absent, $\gamma_t$ outputs collapse to 0),
`cosmological_parameters/{omega_M, omega_nu, h0}`; plus the
single-kernel NFW lookup tables under `data/nfw_off_center/` loaded once
at construction ($c = 4$).

| knob | production value | meaning |
|---|---|---|
| `lambda_bin`, `zo_low`, `zo_high`, `radii` | 120-entry wall (12 bins × 10 $R$) | zipped grid axes |
| `zt_low`, `zt_high` | 0.10, 0.75 | true-$z$ limits |
| `lnm_low`, `lnm_high` | 29.9336, 35.6814 | $\ln M$ GL limits |
| `n_lnm` | 24 | GL nodes in $\ln M$ |
| `n_per_seg` | 10 (default 30) | log-GL nodes per $\theta$ segment |
| `n_zring`, `n_zouter` | 20, 20 | ring / wing $z$ nodes (as in `b_sel_marg`) |
| `n_zt_ref` | 80 | reference-$z$ cache resolution |
| `R_max_cMpch` | 30.0 | sets $\theta_{\max} = R_{\max}/D_A$ |
| `lob_centers` | 25 37.5 52.5 130 (default) | richness centres for $R_{\rm excl}$, $\theta_\lambda$ |
| `theta_breakpoints` | (unset) | optional extra breakpoints, radians |

### Outputs

- `shear_prj/{vals, rnd, cl}` — each length 120: total, the "1"
  (random-point) piece, and the $b\,b_{\rm sel}\,\xi_{\rm NL}$
  (clustering) piece.
- Sibling wrappers around the same core publish
  `sigma_prj/{vals, rnd, cl}` (`SigmaPrjEvaluator`) and
  `dsigma_prj/{vals, rnd, cl}` (`DSigmaPrjEvaluator`) for diagnostics.

### Numerical recipe

Per unique $(\lambda^{\rm ob}, z^{\rm ob})$ slice (grouping all $R$ on
that slice): build one log-GL $\theta$ grid on segments split at the
sorted, deduped breakpoint set
$\{\theta_{\rm lo}, \theta_{\rm excl,o}, \theta_R(R_k) = R_k/D_A,
\theta_\lambda, 2\theta_\lambda, \theta_{\max}\}$ with `n_per_seg` nodes
each — the per-$R$ breakpoint lands nodes on each $\Sigma_{\rm mis}$
peak and cuts the residual $\gtrsim 10\times$ vs a single panel of equal
node count. Inner $z$ uses the same ring + fg/bg $\log|\Delta\chi|$ grid
as `b_sel_marg`; inner $\ln M$ is fixed GL. Exclusion is the per-$z$ LoS
slab $\theta > \theta_{\rm excl}(z)$. The $(z, \ln M)$ integrals are
pre-contracted into per-slice `rnd`/`cl` caches so the outer $\theta$
loop is multiply-adds plus the $\Sigma_{\rm mis}/\Delta\Sigma_{\rm mis}$
table lookup; $\Sigma^{\rm prj}$, $\Delta\Sigma^{\rm prj}$ and
$\gamma_t^{\rm prj}$ come out of one pass ($\sim 250$ ms per sample for
the 120-point wall).

### Implementation

- `src/models/sigma_prj_t.hh` — `sp_detail::ShearPrjCore` (shared
  cache/query backbone, `sp_detail::build_theta_grid`) with thin
  wrappers `ShearPrjEvaluator` / `SigmaPrjEvaluator` /
  `DSigmaPrjEvaluator`, each exporting
  `output_names() = {vals, rnd, cl}`;
- `src/modules/sigma_prj_cpu/ShearPrjEvaluator.cc` (plus
  `SigmaPrjEvaluator.cc`, `DSigmaPrjEvaluator.cc`) —
  `DEFINE_COSMOSIS_SCALAR_EVALUATOR_MODULE(ShearPrjEvaluator)`;
- validation backends sharing the same integrand: **ShearPrjGsl**
  (`ShearPrjGsl.cc`) drives the outer $z$ integral with GSL QAGP
  (explicit breakpoint at $z = z^{\rm ob}$ for the
  $\xi_{\rm NL}(\Delta\chi)$ cusp; knobs `eps_rel`, `eps_abs`,
  `max_eval`, `theta_low`, `theta_high`), writing `sigma_prj_gsl` /
  `dsigma_prj_gsl` / `shear_prj_gsl`; **ShearPrjCuhre**
  (`ShearPrjCuhre.cc`) runs adaptive Cuhre on the inner $(z, \ln M)$
  with the same breakpoint $\theta$ grid, writing `*_cuhre` sections,
  $\sim 30\times$ slower — both off by default and safe to run alongside
  the production evaluator for head-to-head diffs.

### Validation

Mirrors `richness_selection/sigma_prj.py` line-for-line; after the
April–May 2026 audit, $\Sigma^{\rm prj}$ and $\Delta\Sigma^{\rm prj}$
agree with Python to $\leq 1.4\%$ on both `rnd` and `cl` across the full
$R \in [0.1, 30]\,h^{-1}$Mpc wall
(`docs/sigma_prj_py_vs_cpp_report.tex`; residuals traced to
NFW-convention differences, not the projection integrals). Convergence
regression against `ShearPrjCuhre` / `ShearPrjGsl`; Py↔C++ harnesses
(`compare_*_py_vs_cpp.py`) live in
`$PSCRATCH/github/RichnessSelection/validations/`.

*Source: `src/models/sigma_prj_t.hh`, `src/modules/sigma_prj_cpu/`,
`docs/pipeline_modules.tex` §red_shear_prj, the reference ini.*

## likelihoods — Gaussian likelihood

### Scientific definition

Gaussian likelihood comparing the fiducial mock data vector against the
two Costanzi-2026 observables — 12 number counts and 120
tangential-shear points — with the shear theory assembled as a linear sum
of the bin-averaged one-halo term and the projection term:

$$
\begin{aligned}
\gamma_t^{\mathrm{theory}}(R \mid i, j)
  &= \frac{\mathtt{shear1hmissel/vals}}{\mathtt{numcountssel/vals}}
   + \mathtt{shear\_prj/vals}, \\
\log L &= -\tfrac12 \sum_{\mathrm{obs} \in \{\mathrm{NC},\,\mathrm{Shear}\}}
  \delta_{\mathrm{obs}}^{\mathsf T}\, C_{\mathrm{obs}}^{-1}\, \delta_{\mathrm{obs}},
\qquad \delta = d - t .
\end{aligned}
$$

The linear addition of the two shear pieces is valid only because the
reduced-shear $1/(1 - \Sigma\,\Sigma_{\mathrm{crit}}^{-1})$ denominator
was dropped from both integrands (2026-05-11). `shear1hmissel/vals` is
$N_i$-weighted ($N_i[\gamma]$), so dividing by the `NumCountsSel`
integral $N_i[1]$ gives the per-cluster average
$\langle\gamma_t^{1h}\rangle_i$; `shear_prj/vals` is already per-bin.

### Inputs

DataBlock reads: `numcountssel/vals` (12,), `shear1hmissel/vals` (120,),
`shear_prj/vals` (120,). Data/covariance from an `.npz` with keys
`data_NC`, `invcov_NC`, `data_Shear`, `invcov_Shear`; each `invcov` may
be 1-D (diagonal, default) or 2-D dense.

| knob | production value | meaning |
|---|---|---|
| `filename` | `mock_dv_cp_camb.npz` (smoke) / `mock_dv_widePlanck_jkcov.npz` (production MCMC) | mock data vector + inverse covariance |
| `log_space` | F (default) | Gaussian in $y=\ln(\mathrm{obs})$ with delta-method invcov $C_y^{-1}[i,j] = d_i d_j C^{-1}[i,j]$ linearised about the fixed data |
| `verbose` | F | per-sample logL breakdown |

### Outputs

- `likelihoods/likelihoods_like` — scalar $\log L$ (name matched by
  `likelihoods = likelihoods` in the `[pipeline]` block).

### Numerical recipe

1. Setup: load and shape-check the data vector (NC = 12; Shear = 12 bins
   × 10 radii = 120, radius fastest within each bin, matching the
   `shear_prj` wall ordering); if `log_space`, precompute $\ln d$ and the
   delta-method inverse covariance once (constant across the chain, so
   $\det C_y$ and the Jacobian drop out of MCMC).
2. Execute: NC χ² directly on `numcountssel/vals`; shear theory via
   `np.repeat(NC, 10)` tiling, zero-guarded division, plus
   `shear_prj/vals`; χ² as `sum(d²·ic)` (diagonal) or `d @ ic @ d`
   (dense). In log space, theory is floored at $10^{-300}$ so a
   transient non-positive prediction yields a large finite χ² instead of
   a crash.

### Implementation

File: `y3_buzzard/likelihood_cp.py`. Key functions: `_shear_theory`,
`_chi2`, `_to_log_space`, `_residual`, `setup`, `execute`. Layout
constants: `_NC_N_BINS = 12`, `_SHEAR_N_R = 10`, `_SHEAR_N = 120`.

### Validation

Hard size assertions on every datablock read and on the loaded `.npz`
(mismatched vector or covariance shapes abort at setup, not mid-chain).
Closure property: at theory ≡ data, χ² = 0 in both linear and log space,
so fiducial-recovery tests are unaffected by the `log_space` switch; to
first order around fiducial the two chi-squareds agree, only the
off-fiducial posterior geometry changes.

*Source: `y3_buzzard/likelihood_cp.py`, the reference ini,
`docs/pipeline_modules.tex` §observables.*

## Historical appendix

```{toctree}
:maxdepth: 1

historical
```
