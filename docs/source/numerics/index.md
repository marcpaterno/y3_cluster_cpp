# Numerical recipes

The central technical part of this documentation. For every observable the
pattern is: the **full mathematical integral first**, then how the
production implementation reduces it. Where a further factorisation is
planned, the documentation states the **expected factorisation first**
(the spec); the code is then brought to match it, and deltas are flagged.

Primary sources: `docs/pipeline_modules.tex` (algorithms, timing audit,
quadrature knob cheat-sheet), `docs/shear1h_radial_factorization.tex`, and
`RichnessSelection/docs/sigma_prj_refactor.md`.

## General integral structure

The recurring population integral is, schematically,

$$\int dz \int dM \int d\lambda_{\rm true}\;
\mathcal{K}(M, z, \lambda_{\rm true}; \ldots).$$

For each observable the documentation states which dimensions appear in
the full physical expression; which integrations are performed
**analytically** (closed-form richness kernels $\mathcal{K}_i$, $K_j$; the
sigmoid closure in `bsel`); which are **precomputed** per sample
($S_{ij}$ tables, miscentering convolutions); which are **interpolated**;
and which remain **runtime numerical integrations** (typically only the
mass integral).

## Fixed Gauss–Legendre integration (production default)

The Costanzi-2026 branch of the pipeline (`sel_function`, `b_sel_marg`,
`bsel`, `shear_prj`) evaluates its integrals on **fixed Gauss–Legendre
(GL) grids** rather than with adaptive quadrature. A GL rule with $N$
nodes $x_k$ and weights $w_k$ on $[-1, 1]$ is mapped onto each physical
interval $[a, b]$ by the affine transform

$$
t_k = \tfrac{1}{2}(b-a)\,x_k + \tfrac{1}{2}(b+a),
\qquad
\int_a^b f(t)\,dt \approx \tfrac{1}{2}(b-a) \sum_k w_k\, f(t_k).
$$

Two variants of the transform are used:

- **Linear GL** on a fixed bracket — e.g. the $\ln M$,
  $\lambda^{\rm tr}$, and $z$ axes.
- **Log-GL** — GL applied in $\ln\theta$ (nodes placed logarithmically),
  used for the $\theta$ axis of `shear_prj`, where the integrand spans
  decades in angle. The $z$ axis of `shear_prj` and `b_sel_marg` uses a
  "ring + foreground/background" split with log-$|\Delta\chi|$ GL
  spacing: `n_zring` nodes concentrated near $z^{\rm ob}$ and `n_zouter`
  nodes covering the fore/background.

Where GL is used, per module and axis (production node counts from the
reference ini):

| Module | Axis | Rule | Nodes |
|---|---|---|---|
| `sel_function` | $\lambda^{\rm tr}$ | GL on adaptive bracket $[\mu_{\rm eff} - L_\lambda \sigma_{\rm eff},\, \mu_{\rm eff} + L_\lambda \sigma_{\rm eff}]$ per $(\ln M, z)$ | $N_q = 32$ |
| `NumCountsSel` / `Shear1hMisSel` | $z$ | fixed GL, contracted into $W_{ij}(\ln M)$ once per sample | 64 |
| `NumCountsSel` / `Shear1hMisSel` | $\ln M$ | fixed GL (one 1-D sum per wall point) | 96 |
| `b_sel_marg` | $\theta$ | GL, interval split at $\theta_\lambda$ (`n_theta`/2 nodes per segment) | 10 |
| `b_sel_marg` | $z$ | ring + fg/bg GL | 20 + 20 |
| `b_sel_marg` | $\ln M$ | fixed GL (pre-contracted into a $(N_{\rm lt}, N_z)$ cache) | 24 |
| `b_sel_marg` | $\lambda^{\rm tr}$ | fixed GL | 60 |
| `bsel` | $\lambda^{\rm tr}$ | GL on $[1,\, 3\,\lambda^{\rm ob}]$ | 128 |
| `shear_prj` | $\theta$ | log-GL per breakpoint segment | 10 per segment |
| `shear_prj` | $z$ | ring + fg/bg log-$|\Delta\chi|$ GL | 20 + 20 |
| `shear_prj` | $\ln M$ | fixed GL | 24 |

Panels are **hand-placed at integrand features** rather than refined
adaptively. `b_sel_marg` splits its $\theta$ interval at $\theta_\lambda$,
which makes `n_theta` $= 10$ sufficient. `shear_prj` builds one
$\theta$-grid per $(\lambda^{\rm ob}, z^{\rm ob})$ slice, split at the
sorted breakpoint set

$$
\{\theta_{\rm lo},\ \theta_{\rm excl,o},\ \theta_R(R_k)\ \forall R_k,\
\theta_\lambda,\ 2\theta_\lambda,\ \theta_{\max}\},
$$

covering the three integrand features: the exclusion step at
$\theta_{\rm excl,o}$, the sigmoid transition near
$\theta_0 = \theta_\lambda / 2$, and the $\Sigma_{\rm mis}$ peak at
$\theta_R(R) = R / D_A$. A single GL panel over
$[\theta_{\rm lo}, 2\theta_\lambda]$ under-resolves all three; the
breakpoint split reduces the residual by $\gtrsim 10\times$ at
`n_per_seg` $= 30$ compared to a single-panel 120-node grid.

**The $n_{\ln M}$ convergence study.** The `sel_function` mass grid was
swept as a *whole-pipeline* study, because the `S_stack` grid resolution
feeds directly into how hard the downstream Cuhre-based `NumCountsSel`
must refine:

| $n_{\ln M}$ | `sel_function` | `NumCountsSel` | `Shear1h` | total | $\Delta N_c$ vs 256 |
|---:|---:|---:|---:|---:|---|
| 256 | 0.384 s | 0.083 s | 0.097 s | 1.36 s | ref |
| **192** | **0.224 s** | **0.121 s** | 0.085 s | **1.15 s** | $+0.01\%$ |
| 128 | 0.148 s | 0.273 s | 0.080 s | 1.22 s | $+0.03\%$ |
| 96 | 0.113 s | 0.211 s | 0.076 s | 1.13 s | $+0.07\%$ |
| 64 | 0.072 s | 0.207 s | 0.077 s | 1.11 s | $-4.55\%$ — **GL resonance, avoid** |
| 48 | 0.055 s | 0.262 s | 0.062 s | 1.18 s | $+0.33\%$ |

The optimum is $n_{\ln M} = 192$, *not* the single-module optimum: below
192, Cuhre reclaims the `sel_function` savings by refining more; and
$n_{\ln M} = 64$ hits a **pathological GL resonance** producing a $4.5\%$
error on the counts — note the non-monotone error (48 is fine, 64 is
not), which is why sweeps must bracket the target, not just sample it.
Accuracy at 192 vs 256 is $< 0.05\%$ on `NumCountsSel` and `Shear1hSel`.

**Rationale.** A fixed-GL evaluator performs an input-independent number
of integrand evaluations per sample: its cost is deterministic across the
parameter space. Adaptive integrators (Cuhre) refine to a tolerance, so
their cost depends on the integrand at each MCMC point — exactly the
coupling exposed by the sweep above, where the *same* Cuhre settings cost
0.083 s or 0.273 s depending on upstream grid resolution. For MCMC,
deterministic per-sample cost removes the adaptive-refinement tail
latency. It also permits structural optimisations that adaptivity
forbids: $\theta$-outer loop nests with pre-contracted
$(\lambda^{\rm tr}, z)$ mass caches (`b_sel_marg`), per-slice
$z$/$\ln M$-contracted weights (`shear_prj`, a $\sim 2\times$ win over
the naive $(z, \theta, M)$ nest), and emitting $J = I_2 - I_1$ directly
from its non-negative integrand
$(1 - \sigma(\theta))\,\xi_{\rm NL}\, b\, f_A$ instead of differencing
two nearly-equal GL sums (which loses $\sim 1$ digit to cancellation as
$\sigma(\theta) \to 1$).

*Source: `docs/pipeline_modules.tex` §sel_function, §b_sel_marg,
§red_shear_prj, §Timing & precision optimisation audit (node counts from
the reference ini).*

## Adaptive integration (reference and validation backends)

Three adaptive backends remain in the pipeline, in distinct roles:

**CUBA/Cuhre (CPU, retired for the $N_i[f]$ family, kept for the
brute-force reference).** Until 2026-08 `NumCountsSel` and the one-halo
shear module ran the CUBA Cuhre cubature with `eps_rel`
$= 1.5\times10^{-3}$, `eps_abs` $= 10^{-12}$, `max_eval` $= 10^6$, and
`use_cartesian_product = T`, at 106 ms and 131 ms per sample (mean; the
adaptive tail reached 0.98 s and 4.0 s respectively over $\sim 10^6$ MCMC
realisations). Both now use the fixed-GL evaluator
(`n_operator_sel_gl_t.hh`) with deterministic cost — the Cuhre knobs
still present in the ini sections are **ignored**. Cuhre remains the
driver of the brute-force triple-integral reference
`NumCountsFullScalarIntegrand` (`max_eval` $= 10^7$), which mirrors the
emulator-validation report's configuration so cp_camb-vs-CAMB residuals
stay comparable. The tolerance-driven cost coupling documented in the
$n_{\ln M}$ sweep above is the historical record of why the Cuhre path
was retired.

**GSL QAGP.** `ShearPrjGsl` runs adaptive GSL QAGP over the $z$ axis of
the projection integrand for tight-tolerance diagnostics. Off by default.

**Cuhre variant of the projection.** `ShearPrjCuhre` runs adaptive Cuhre
on the inner $(z, \ln M)$ integral of the same integrand implementation
as the production evaluator. Off by default, $\sim 30\times$ slower than
the fixed-GL `ShearPrjEvaluator` ($\sim 250$ ms/sample); used for
convergence regression tests. All three projection evaluators share one
integrand implementation, so backend comparisons isolate quadrature
error.

**PAGANI (GPU).** The repo also builds PAGANI-backed GPU variants of the
integrands. They are reference/benchmark only, not a production path: on
the `b_sel` wall grid the fixed-GL evaluator takes 0.17 s where the
PAGANI variant takes 208 s — roughly $10^3\times$ slower for these
low-dimensional, feature-dominated integrands.

**Timing-audit trajectory.** The April–May 2026 audit reduced the
per-sample pipeline cost by $\sim 560$ ms ($-33\%$):

| Change | Total per sample |
|---|---:|
| Baseline | **1.71 s** |
| `compute_lensing_2h = F` | 1.31 s |
| `_plob_params` 1-D $z$ fast path | 1.25 s |
| `_cdf_lob` in-place | 1.22 s |
| `n_lnm = 192` | **1.15 s** |

Re-measured on `mock_mcmc_cp_camb.ini` with `Shear1hMisSel` as the
1-halo branch, the total is $\sim 1.14$ s (login node, single thread),
dominated by `shear_prj` (231 ms), `sel_function` (215 ms), `mf_tinker`
(175 ms), and `halo_model` (158 ms).

*Source: `docs/pipeline_modules.tex` §Pipeline topology,
§red_shear_prj ("Three integration backends, same recipe"), §Timing &
precision optimisation audit (Cuhre settings from the reference ini).*

## Production quadrature settings ("the wall of numbers")

All values below are the live settings in the reference production ini.
Mass bounds are given as $\ln M$ with $M$ in $h^{-1} M_\odot$:
$29.9336 = \ln 10^{13}$, $35.6814 = \ln 10^{15.5}$,
$36.7300 \approx \ln 10^{15.95}$, $36.8414 = \ln 10^{16}$.

| Knob | Module | Value | What it controls | Why that value |
|---|---|---|---|---|
| `lam_min` / `lam_max` | `sel_function` | $\{20, 30, 45, 60, 200\}$ edges | richness-bin edges of $K_i$ | DES Y3 richness binning; 4 bins $\times$ 3 $z$ bins = 12-cell wall |
| `zob_min` / `zob_max` | `sel_function` | $\{0.20, 0.35, 0.50, 0.65\}$ edges | photo-$z$ bin edges of $K_j$ | DES Y3 photo-$z$ binning |
| `sigma_z` | `sel_function` | 0.03 | photo-$z$ scatter in the Gaussian kernel $K_j(z)$ | survey photo-$z$ scatter model |
| `zt_low` / `zt_high` | `sel_function` | 0.05 / 0.80 | true-$z$ tabulation range of $S_{ij}$ | brackets the $z^{\rm ob}$ bins $[0.20, 0.65]$ with room for the $\sigma_z = 0.03$ Gaussian tails |
| `lnm_low` / `lnm_high` | `sel_function` | 29.9336 / 36.8414 | $\ln M$ tabulation range | $[10^{13}, 10^{16}]\,h^{-1}M_\odot$ — full cluster mass range; upper edge above every consumer's integration ceiling |
| `n_lnm` | `sel_function` | **192** | $\ln M$ grid size of `S_stack` | whole-pipeline optimum of the 2026-05-07 sweep; $<0.05\%$ shift vs 256; 64 is a GL-resonance pathology (avoid) |
| `n_z` | `sel_function` | 20 | module-internal $z$ grid | ini default (the tabulation grid is `n_z_shared`) |
| `n_z_shared` | `sel_function` | 64 | shared $z$ grid size of `S_stack` | matches the 1-D $z$ fast path of the 8 EMG coefficient splines |
| `L_lam` | `sel_function` | 6.0 | half-width of the $\lambda^{\rm tr}$ GL bracket $[\mu_{\rm eff} \pm L_\lambda \sigma_{\rm eff}]$ | $\pm 6\sigma$ captures the shifted-Poisson $P_{\rm HOD}$ tails |
| `L_z` | `sel_function` | 6.0 | analogous bracket half-width for the $z$ axis | $\pm 6\sigma$ coverage of the Gaussian $K_j$ kernel |
| `N_q` | `sel_function` | 32 | GL nodes in $\lambda^{\rm tr}$ per $(\ln M, z)$ | GL bracket nodes in $\lambda^{\rm tr}$ (knob cheat-sheet) |
| `algorithm`, `eps_rel`, `eps_abs`, `max_eval`, `use_cartesian_product` | `NumCountsSel`, one-halo shear | (legacy) | Cuhre knobs from the retired adaptive path | **ignored** by the fixed-GL evaluators; still honoured by the `NumCountsFullScalarIntegrand` brute-force reference (`eps_rel` $=1.5\times10^{-3}$, `max_eval` $=10^7$) |
| `n_lnm` / `n_z` | `NumCountsSel`, one-halo shear | 96 / 64 (defaults) | fixed-GL node counts | z-contracted weight + 1-D mass sum per wall point |
| `zt_low` / `zt_high` | `NumCountsSel`, one-halo shear | 0.05 / 0.80 | true-$z$ integration range | identical to the `sel_function` tabulation range — no clamped edge queries |
| `lnm_low` / `lnm_high` | `NumCountsSel`, one-halo shear | 29.9336 / 36.7300 | $\ln M$ integration range | $[10^{13}, \sim 10^{15.95}]$; stays strictly inside the `S_stack` mass table ($\ln M \le 36.8414$) |
| `r_perp` | one-halo shear | 10 log-spaced, 0.10–10.0 | radii of the 1-halo shear wall | $R$ grid of the shear datavector, $h^{-1}\mathrm{Mpc}$ |
| `lnm_low` / `lnm_high` | `b_sel_marg` | 29.9336 / 35.6814 | $\ln M$ range of the $P[X]$ operators | $[10^{13}, 10^{15.5}]$ — matches the Python `bsel` reference mass grid |
| `n_lt` | `b_sel_marg` | 60 | $\lambda^{\rm tr}$ GL nodes | already converged at 16 (cheat-sheet); 60 is comfortable margin |
| `n_lnm` | `b_sel_marg` | 24 | $\ln M$ GL nodes (pre-contracted mass cache) | inner mass integral is smooth; contraction makes it cheap |
| `n_theta` | `b_sel_marg` | 10 | $\theta$ GL nodes, split at $\theta_\lambda$ | sufficient *with* the split at the sigmoid feature |
| `n_zring` / `n_zouter` | `b_sel_marg` | 20 / 20 | $z$ nodes near $z^{\rm ob}$ / fore-background | ring resolves the $\xi_{\rm NL}(|\Delta\chi|)$ peak at $z \approx z^{\rm ob}$ |
| `lob` / `zob` | `bsel` | $\{25, 37.5, 52.5, 130\}$ / $\{0.275, 0.425, 0.575\}$ | wall centres | DES Y3 arithmetic bin centres (130 replaces the legacy 100) |
| `n_theta`, `theta_lo`/`theta_hi` | `bsel` | 32 on $[10^{-4}, 5\times10^{-3}]$ rad | $\theta$ grid of the published $b_{\rm sel}(\theta)$ table | diagnostic table only; `shear_prj` evaluates the sigmoid analytically |
| `n_ltr` | `bsel` | 128 | $\lambda^{\rm tr}$ GL nodes of the $(b_{\rm zero}, b_\infty)$ marginalisation | matches the Python `SelBias._marginalised_plateaus` reference |
| `ltr_lo` / `ltr_hi_factor` | `bsel` | 1.0 / 3.0 | marginalisation range $[1,\, 3\,\lambda^{\rm ob}]$ per bin | matches the Python reference range |
| `zt_low` / `zt_high` | `shear_prj` | 0.10 / 0.75 | line-of-sight $z$ integration range | brackets the lens bins; the parabolic $w_z(z, z^{\rm ob})$ kernel gates the interior |
| `lnm_low` / `lnm_high` | `shear_prj` | 29.9336 / 35.6814 | $\ln M$ integration range | $[10^{13}, 10^{15.5}]$, same as `b_sel_marg` — the two Costanzi-2026 operators share a mass convention |
| `R_max_cMpch` | `shear_prj` | 30.0 | upper $\theta$ boundary of the projection integral | comoving cap; validation covers $R \in [0.1, 30]\,h^{-1}\mathrm{Mpc}$ |
| `n_lnm` | `shear_prj` | 24 | $\ln M$ GL nodes (pre-contracted per slice) | smooth inner integral, contracted into per-slice caches |
| `n_per_seg` | `shear_prj` | 10 | log-GL nodes per breakpoint $\theta$ segment | ini production value (module default 30); the feature-aligned breakpoints do the resolving, giving $\gtrsim 10\times$ residual reduction vs a single 120-node panel |
| `n_zring` / `n_zouter` | `shear_prj` | 20 / 20 | ring + fg/bg log-$|\Delta\chi|$ $z$ nodes | same ring structure as `b_sel_marg` |
| `n_zt_ref` | `shear_prj` | 80 | reference $z_t$ grid size of the evaluator | ini production value |
| `lambda_bin` / `zo_*` / `radii` | `shear_prj` | 120-point wall | $(\lambda^{\rm ob}, z^{\rm ob}, R)$ grid | 12 bins $\times$ 10 radii, $R \in [0.10, 10.0]$ matching the one-halo `r_perp` grid |

*Source: `docs/pipeline_modules.tex` §sel_function, §b_sel_marg, §bsel,
§red_shear_prj, §Knob cheat-sheet (values from the reference ini).*

## Interpolation and precomputed tables

The C++ evaluators never compute physics primitives in their inner loops;
everything slowly varying is tabulated once per sample (or once per
process) and served through interpolators.

**Per-sample tables published on the datablock:**

- $S_{ij}(\ln M, z)$ — `sel_function/S_stack`, shape
  $(N_{\rm bin}, N_z, N_{\ln M}) = (12, 64, 192)$ on linear grids
  $\ln M \in [29.9336, 36.8414]$, $z \in [0.05, 0.80]$. Served to
  `NumCountsSel` and the 1-halo shear module via `Interp2D`. Grid size is
  the whole-pipeline optimum (see the $n_{\ln M}$ study): the
  interpolation error of the tabulated selection function feeds directly
  into how hard the downstream adaptive integrator refines.
- $b(M, z)$ — `haloModel/bias`, Tinker-2010 bias on a 100-node log mass
  grid $[10^{12}, 10^{16}]$ with $z$ from the power-spectrum grid
  ($n_z = 50$); consumed directly via `Interp2D` by `b_sel_marg` and
  `shear_prj`. (There is no compiled halo-bias type; the table *is* the
  interface.)
- $\xi_{\rm NL}(r, z)$ — `xi_nl`, tabulated over
  $r \in [10^{-3}, 10^3]\,h^{-1}\mathrm{Mpc}$ on 128 nodes per $z$ slice;
  queried at $|\Delta\chi(\theta, z)|$ in the $\theta$-outer loops.
- Centred NFW lensing — `haloModel/Sigma_nfw`, `dSigma_nfw` on
  $R_\perp \in [0.05, 10.0]$, 128 bins (Child-18 concentration); the
  `Shear1hMisSel` weight reads it as `dSigma_nfw.clamp(R, lnM)`.
- $(B_{\rm small}, B_{\rm large})(z^{\rm ob})$ — `shear_prj` **linearly
  interpolates** the `b_sel_marginalised` scalars to the exact
  $z^{\rm ob}$ of each slice, then evaluates
  $b_{\rm sel}(\theta) = B_{\rm small}(1 - \sigma(\theta)) +
  B_{\rm large}\,\sigma(\theta)$ **analytically**. It deliberately does
  *not* spline the legacy 32-node $b_{\rm sel}(\theta)$ table: no spline
  ringing, and no nearest-$z^{\rm ob}$ tie-break bug.
- HMF — `HMF_t` wraps the `mass_function/dndlnmh` spline. The
  constructor rescales the mass axis internally (CosmoSIS `m_h` is
  stored in $\Omega_m \cdot h^{-1} M_\odot$ Komatsu-CRL units), so every
  consumer queries with raw $\ln M_{\rm phys}$ and must set
  `ln_mass_shift` $= 0$ — a second shift is a double-count that moves
  the mass axis by $\approx 0.6$ dex.

**Process-lifetime disk tables (miscentered NFW):**

The azimuthal-plus-radial miscentering convolution

$$
\begin{aligned}
\Delta\Sigma_{\rm mis}(R, M; \tau_{\rm mis})
  &= \int_0^\infty dR_{\rm mis}\, P(R_{\rm mis} \,|\, \tau_{\rm mis}, R_\lambda)\\
  &\quad\times \frac{1}{2\pi}\int_0^{2\pi} d\varphi\,
     \Delta\Sigma_{\rm NFW}\!\big(\sqrt{R^2 + R_{\rm mis}^2 + 2 R R_{\rm mis}\cos\varphi}\,\big)
\end{aligned}
$$

is far too expensive for inner loops ($\sim 10^5$ integrand
nodes/sample), so both consumers read pre-computed $1000 \times 1000$
**log–log grids** in $(x = R/r_s,\ x_{\rm mis} = R_{\rm mis}/r_s)$ from
`data/nfw_off_center/`, storing
$\log[\Delta\Sigma_{\rm mis} / (2\, r_s\, \delta_c\, \rho_{\rm crit})]$
for the `single` ($\delta$-function offset; used by `shear_prj`, where
the offset *is* the $\theta$ variable) and `gamma` (Rayleigh-shaped
kernel; used by the 1-halo mixture) kernels. Tables load once at module
construction; at evaluation time $r_s = r_{200}/c$ (with $c = 4$)
rescales the dimensionless profile — no on-the-fly convolution. The
per-node `Interp2D::clamp` lookup costs $\sim 32$ ms/sample relative to
the centred-only branch (131 ms vs 99 ms).

**Extrapolation behavior.** Out-of-range queries are **clamped** to the
table edge (`Interp2D::clamp`, `dSigma_nfw.clamp`) rather than
extrapolated. The production ini keeps integration ranges strictly inside
tabulation ranges (e.g. $\ln M \le 36.7300$ vs a table ceiling of
$36.8414$; identical $z_t$ ranges for `sel_function` and its consumers)
so clamping is a safety net, not an active code path.

**Accuracy.** With these tables and the production grids, the C++
pipeline agrees with the independent Python reference to $\le 1.4\%$ on
$\Sigma^{\rm prj}, \Delta\Sigma^{\rm prj}$ (both `rnd` and `cl` pieces)
across $R \in [0.1, 30]\,h^{-1}\mathrm{Mpc}$, $\le 1\%$ on
$(P_1, I_1, J)$, and $\sim 1\%$ on $(B_{\rm small}, B_{\rm large})$.

*Source: `docs/pipeline_modules.tex` §Shared conventions, §halo_model,
§sel_function, §red_shear_prj, §Miscentering selection on
$\Delta\Sigma$, §Timing & precision optimisation audit (grid ranges from
the reference ini).*

## Unit and convention traps

Hard-won conventions that every new integrand must respect:

- **chi vs $R_\lambda$ units**: `distances/d_c` returns $\chi(z)$ in
  **Mpc**; $R_\lambda = (\lambda/100)^{0.2}$ is in **cMpc/h**. Consumers
  must multiply $\chi$ by $h_0$ to put both on the same axis.
- **HMF mass axis**: the HMF table stores $dn/d\ln M$ against
  $\ln\!\left(m_h(\Omega_m - \Omega_\nu)\right)$. Because the loader
  already applies the rescaling, `ln_mass_shift_` must be **0.0** in all
  evaluators; a raw $\ln M$ query returns the value at the wrong mass
  (≈ 0.6 dex shift) and silently produces wrong HMF values.
- **$\lambda_{\rm true}$ range**: Costanzi-2026 operators integrate
  $\lambda_{\rm true}$ on $(0, \lambda_{\rm ob,centre}]$ *per grid
  point*, not on a fixed global interval.
- **$B_i$ and $\Omega(z)$ placement**: both are deliberately absent from
  the $P[X]$ operator — they cancel in every downstream $b_{\rm sel}$
  ratio; $\Omega(z)$ is likewise excluded from surface-density
  observables (see {doc}`../systematics/index`).

## Implementation mapping

The equation → algorithm → code mapping for every production module —
DataBlock inputs/outputs, ini knobs, source files, and classes — lives in
the {doc}`module reference <../modules/index>`. The GPU (PAGANI) variants
are reference benchmarks only (see the adaptive-integration section
above) and the historical GPU suites are catalogued in the
{doc}`status appendix <../modules/historical>`.
