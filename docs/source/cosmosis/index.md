# The CosmoSIS pipeline

How `y3_cluster_cpp` fits into a complete CosmoSIS pipeline, and the full
trace of the reference configuration.

## Three categories of modules

1. **CosmoSIS Standard Library** — parameter consistency
   (`consistency`), linear growth (`GrowthFactor`), the Tinker mass
   function (`MfTinker`), distances. The full-fidelity $P(k)$ path is
   CAMB; the reference pipeline replaces it with the `cp_camb`
   CosmoPower emulator.
2. **This analysis's own modules** — the C++ `.so` files under
   `release-build/src/modules/…` and the Python modules under
   `y3_buzzard/` and `src/modules/…/*.py`, documented in
   {doc}`../modules/index`.
3. **External dependencies** — CUBA/cubacpp, PAGANI
   (`gpuintegration`), GSL, `cluster_toolkit`, and the CUDA toolkit for
   the GPU variants.

## The reference configuration

`des-cluster-nersc/cosmosis-models/mock_mcmc_cp_camb.ini` @ `5055375`
(2026-08-09). Module list:

```text
consistency  GrowthFactor  cp_camb  MfTinker  halo_model
average_sigma_crit_inv  sel_function
NumCountsSel  Shear1hMisSel
b_sel_marg  bsel
shear_prj
likelihoods
```

Key facts:

- **Samplers**: `test` for a single smoke sample; `apriori` for prior
  coverage; `emcee` (64 walkers) and `polychord` (500 live points) for
  production. Closure requirement: $\log L \approx 0$ at the fiducial
  point, and the posterior recovers the fiducial HOD parameters.
- **Values**: `mock_mcmc_widePlanck_values.ini` — 5 cosmology + 5 HOD
  parameters varied.
- **Data vector**: `mock_dv_widePlanck_jkcov.npz` — 12 number counts +
  120 shear points; Poisson number-count covariance and the Buzzard
  jackknife shear covariance.
- **Naming convention**: C++ modules use CamelCase section names, Python
  modules snake_case — except modules whose source hardcodes
  `module_label()` (`b_sel_marg`, `shear_prj`), which keep their legacy
  lowercase names.
- The copy of this file inside `y3_cluster_cpp/cosmosis-models/` is the
  **single-sample smoke variant** (test sampler, centered-only
  `Shear1hSel`).

Run pipelines from the sibling repositories with `Y3_CLUSTER_CPP_DIR`
pointing at this tree, so the built `.so` files under
`release-build/src/modules/…` resolve. Run management (sbatch scripts,
chain handling, mock-data-vector generation) lives in
`des-cluster-nersc`; Python↔C++ validation harnesses in
`RichnessSelection`; the $P(k)$ emulator training in `camb-emulator`.

## DataBlock trace

`docs/figs/real_pipeline_extract.ini` is a trimmed copy of the reference
configuration that runs the first nine stages (consistency →
`Shear1hMisSel`) under the test sampler, dumping every DataBlock section
to `real_pipeline_extract_output/`. Shapes below are read directly from
those dumps; stages past `Shear1hMisSel` are described from the module
sources since the extract stops before them.

| Stage (module) | Sections written | Key contents (names + shapes) |
|---|---|---|
| `consistency` (CSL) | `cosmological_parameters` (completed) | Derived-parameter closure: `h0`/`hubble`, `omega_m`, `omega_b`, `omega_c`, `omega_nu`, `omega_k`, `omega_lambda`, `ommh2`, `ombh2`, `omch2`, `sigma8`, `n_s`, `mnu`, `w`, `wa`, `baryon_fraction` (scalars) |
| `GrowthFactor` (CSL) | `growth_parameters` | `z`, `a` (406 pts, $z \in [0, 4.05]$, `dz` = 0.01), `d_z` (linear growth), `f_z` (growth rate) |
| `cp_camb` | `matter_power_lin`, `cdm_baryon_power_lin`, `distances` | `z` (50, $z \in [0,4]$), `k_h` (506), `p_k` (50 × 506); same grid for the no-neutrino $P_{cb}$; `distances/{z, a, d_a, d_c, d_l, d_m, h, mu, nz}` (50 pts) |
| `MfTinker` (CSL) | `mass_function` | `m_h` (969), `z` (50), `dndlnmh` (50 × 969), plus `r_h`, `dndlnrh`; mass axis is $M_{\rm phys}/\Omega_m$ (the `HMF_t` consumers apply the $\ln(\Omega_m-\Omega_\nu)$ shift internally) |
| `halo_model` | `halomodel`, `xi_nl` | `halomodel/{lnm (100), m_h (100), z (50), bias (50 × 100), concentration (100), rhoc, r_sigma (128), k (506), sigma_nfw (100 × 128), dsigma_nfw (100 × 128), hubble_shift, scale_shift}`; `xi_nl/{r (128), z (50), xi_nl (50 × 128)}`. With `compute_lensing_2h = F` no `Sigma_hh`/`dSigma_hh`/`Wp_hh` keys appear. `bias` is the sole halo-bias source downstream |
| `prj_params` (smoke ini only) | `plob_ltr_params` | `z` (15 nodes, $z \in [0.10, 0.80]$) + the 8 EMG coefficient arrays `a_tau, b_tau, a_mu, b_mu, a_sig, b_sig, a_fprj, b_fprj` (each 15). The production configuration drops this module — `sel_function` and `bsel` fall back to the identical table embedded in `y3_buzzard/prj_params.py` |
| `average_sigma_crit_inv` | `average_sigma_crit_inv` | `zlense` (50, $z \in [0.05, 0.80]$), `sci_average` (50). Under `unity = T` the dumped `sci_average` is identically 1.0; otherwise it integrates the `test_cluster_Y1.fits` source $n(z)$ against `distances/d_a` |
| `sel_function` | `sel_function` | `lnm` (192), `z` (64), `s_stack` with shape (12, 64, 192) — the per-bin selection surface $S_{ij}(\ln M, z)$ |
| `NumCountsSel` | `numcountssel` | `vals` (12): expected counts $N_i[1]$ per bin |
| `Shear1hMisSel` | `shear1hmissel` | `vals` (120 = 12 bins × 10 `r_perp`): miscentering-weighted 1-halo stack numerators $N_i[\gamma_t^{1h,\rm full}](R)$; richness bin resolved as `bin_index % 4` |
| `b_sel_marg` | `b_sel_marg_P1`, `b_sel_marg_I1`, `b_sel_marg_J` | Each: `vals` on the 12-point $(z_{\rm ob}, \lambda)$ wall + grid echo `zo_low`, `zo_high`, `lambda_bin`; $J = I_2 - I_1$ is computed directly to avoid cancellation at large $\theta$ |
| `bsel` | `b_sel_marginalised` | `lob` (4), `zob` (3), `theta` (32), `vals` (4 × 3 × 32) — the $\lambda_{\rm true}$-marginalised $b_{\rm sel}(\theta)$ lookup — plus `b_eff`, `b_small`, `b_large` |
| `shear_prj` | `shear_prj` | `vals`, `rnd`, `cl` (each 120): projection shear split as total = rnd + cl (uncorrelated random-point term + clustered $\xi_{\rm NL}$-weighted, $b_{\rm sel}$-modulated term) |
| `likelihoods` | `likelihoods` | `likelihoods_like` (scalar Gaussian $\log L$); assembles the theory vector from `numcountssel/vals`, `shear1hmissel/vals`, `shear_prj/vals` against the mock `.npz` |

Two notes on the dump set. First, the extract output also contains
`cluster_abundance` (`hmf_q`, `hmf_s`), `cluster_mor` (`alpha`,
`epsilon`, `log10_mmin`, `log10_ratio`, `sigma_lambda`), and `photoz`
(`delta_z`) — sampler-parameter sections seeded from the values file,
not module outputs; they are read by `sel_function`, `bsel`, and the
MOR/HMF models. Second, every dumped section remains live in the current
production pipeline (none are orphaned), though
`average_sigma_crit_inv/sci_average` is unity-valued only under
`unity = T`.
