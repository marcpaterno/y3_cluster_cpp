# `y3_cluster_cpp` code reference

A source-tree reference for `src/models/` (header-only physics terms) and
`src/modules/` (CosmoSIS `.so` modules). Inventory-style: one row per file,
grouped by physics domain (models) or build status (modules).

This is **not** the pipeline reference. For wired-pipeline algorithms,
datablock contracts, and Costanzi-2026 numerical recipes, see
[`pipeline_modules.tex`](pipeline_modules.tex) (PDF: `pipeline_modules.pdf`).
For build instructions, see [`../BUILDING.md`](../BUILDING.md).

Conventions:

- `**module**` — a CosmoSIS `.so` produced by a directory under `src/modules/`.
- `` `section/key` `` — a CosmoSIS datablock section/key.
- `` `path/to/file` `` — a path under the repo root.
- "GPU?" column = `✓` if a matching `.cuh` device variant exists.

---

## Contents

1. [Overview](#1-overview)
2. [Most-updated modules (Costanzi-2026 production path)](#2-most-updated-modules-costanzi-2026-production-path)
3. [`src/models` — physics terms](#3-srcmodels--physics-terms)
4. [`src/modules` — CosmoSIS `.so` registry](#4-srcmodules--cosmosis-so-registry)
5. [Module-macro primer](#5-module-macro-primer)
6. [Appendix A — Retired / not built / unused](#appendix-a--retired--not-built--unused)
7. [Appendix B — Python-only modules](#appendix-b--python-only-modules-under-srcmodules)

---

## 1. Overview

The C++/CUDA code in this repo is organised around two ideas:

- **`src/models/`** holds *integrand terms* — header-only template structs
  that evaluate one physical quantity at a point (a MOR probability, an
  HMF amplitude, an NFW Σ, a richness kernel, …). Each constructor takes
  a `cosmosis::DataBlock&` and reads its inputs from named sections.
- **`src/modules/`** holds *CosmoSIS modules* — thin `.cc` drivers that
  instantiate an integrand template, pick an integrator (CUBA / Cuhre /
  PAGANI / fixed-GL evaluator), and emit the `setup`/`execute`/`cleanup`
  C-ABI via a `DEFINE_COSMOSIS_*_MODULE` macro from
  `src/utils/module_macros.hh` (see §5). Each subdirectory builds one or
  more `.so` files registered in `src/modules/CMakeLists.txt`.

The current production path is the **Costanzi-2026** pipeline: closed-form
selection-bias `N`-operators (`num_counts_sel`) plus the fixed-GL evaluators
`b_sel_marg_cpu` (P₁, I₁, J co-pass) and `sigma_prj_cpu` (Σ_prj + ΔΣ_prj +
γ_t^prj co-pass). The PAGANI/CUDA variants under
`b_sel_marg_cpu/*PaganiIntegrand.cc`, `sigma_park_y1`, `sigma_buzzard_y3`,
and `gt_*_gpu` are kept as benchmarks or historical Y1 work — see §4.

---

## 2. Most-updated modules (Costanzi-2026 production path)

| Module (`.so`) | Source | Computes | See also |
|---|---|---|---|
| **BSelMargIntegrand** | `src/modules/b_sel_marg_cpu/BSelMargIntegrand.cc` | Fixed-GL evaluator co-computing `b_sel_marg/{P1,I1,J}` on a common `(zob_low, zob_high, λ_bin)` grid. `J` is computed directly to avoid catastrophic cancellation against `I1·I2`. Replaces three PAGANI integrands. | `pipeline_modules.tex` §8 *b_sel_marg* (line 582) |
| **ShearPrjEvaluator** | `src/modules/sigma_prj_cpu/ShearPrjEvaluator.cc` | Fixed-GL evaluator: single outer θ loop with log-GL breakpoints, inner ring + foreground/background `(z, lnM)` quadrature, co-computes `sigma_prj`, `dsigma_prj`, and `shear_prj` outputs. Replaces a three-module split. Compile flags: `-ffast-math -fopenmp-simd -march=native`. | `pipeline_modules.tex` §10 *red_shear_prj* (line 742) |
| **NumCountsSel + 10 weight variants** | `src/modules/num_counts_sel/{NumCounts,MassWeighted,BiasWeighted,Sigma1h,SigmaTot,DSigma1h,DSigma1hMis,DSigmaTot,Shear1h,Shear1hMis,ShearTot}.cc` | Closed-form selection-bias `N`-operators `N_i[f]`. Each `.cc` is a thin (≈10-line) typedef of `NOperatorSelScalar<W>` or `NOperatorSelRadial<W>` from `src/models/n_operator_sel_t.hh`; the weight `W` lives in `lensing_weights.hh`. Consumes `S_ij(lnM, z)` from the `sel_function` Python module. | `pipeline_modules.tex` §3.1 *N_i[f] operator* (line 250); §11 *Miscentering selection* (line 860) |
| **halo_model** (Python sibling) | (driven from sibling repos; `.py` writes datablock `halomodel/bias` and `halomodel/xinl`) | Publishes `b(M, z)` (`Interp2D`) and `ξ_NL(r, z)` consumed by `sigma_prj_cpu` and downstream lensing 2-halo integrals. | `pipeline_modules.tex` §6 *halo_model* (line 418) |
| **sel_function** (Python sibling) | (`.py` only — no `.so`) | Tabulates `S_ij(lnM, z)` consumed by every `*Sel` module. EMG + HOD math, GL bracket nodes, 1-D z fast-path. | `pipeline_modules.tex` §7 *sel_function* (line 486) |

Cross-cutting note: on the b_sel wall grid the fixed-GL evaluator
(BSelMargIntegrand) runs in ≈0.17 s vs ≈208 s for the PAGANI integrands —
roughly 10³× — so the PAGANI `.cc` files in `b_sel_marg_cpu/` are
benchmark-only, not production.

---

## 3. `src/models` — physics terms

All headers are header-only template structs in `namespace y3_cluster`.
"DataBlock" lists the section the constructor reads (or *static* if it uses
a frozen embedded table). "GPU?" = ✓ when a matching `.cuh` exists.

### 3.1 Cosmological basics

| Header | Type | Physics | DataBlock | `operator()` | GPU? |
|---|---|---|---|---|---|
| `ez.hh` | `EZ` | E(z) = H(z)/H₀ | `cosmological_parameters` | `(z) → double` | ✓ |
| `ez_sq.hh` | `EZ_sq` | E(z)² | `cosmological_parameters` | `(z) → double` | ✓ |
| `dv_do_dz_t.hh` | `DV_DO_DZ_t` | dV / dΩ dz comoving volume element | `distances` | `(zt) → double` | ✓ |
| `omega_z_des.hh` | `OMEGA_Z_DES` | Ω(z) DES survey area | `survey` | `(zt) → double` | ✓ |
| `omega_z_sdss.hh` | `OMEGA_Z_SDSS` | Ω(z) SDSS survey area | `survey` | `(zt) → double` | — |
| `sigma_crit_inverse_t.hh` | `sigma_crit_inv` | Σ_crit⁻¹(z_lens, z_src) | `distances` | `(zt, zs) → double` | — |
| `average_sci_t.hh` | `AVERAGE_SCI_t` | Source-redshift-averaged ⟨Σ_crit⁻¹⟩(z_lens) | `pz_source` | `(zt) → double` | — |
| `weighted_sigma_crit_inv.hh` | `WeightedSigmaCritInv` | Weighted ⟨Σ_crit⁻¹⟩ | `pz_source` | `(zt, bin) → double` | — |
| `angle_to_dist_t.hh` | `ANGLE_TO_DIST_t` | θ → physical R(z) | `distances` | `(theta, zt) → double` | — |

### 3.2 Halo model

| Header | Type | Physics | DataBlock | `operator()` | GPU? |
|---|---|---|---|---|---|
| `hmf_t.hh` | `HMF_t` | n(M, z) Tinker mass function | `mass_function` | `(lnM, zt) → double` | ✓ |
| `nfw_sigma_t.hh` | `nfw_sigma` | NFW Σ(R, M) | static | `(r, M) → double` | — |
| `nfw_sigma_mis.hh` | `NFW_SIGMA_MIS` | Off-centre NFW Σ_mis(R, M, R_off) (active) | `nfw_off_center` | `(r, M, R_off) → double` | ✓ |
| `nfw_dsigma_mis.hh` | `NFW_DSIGMA_MIS` | Off-centre NFW ΔΣ_mis (active) | `nfw_off_center` | `(r, M, R_off) → double` | ✓ |
| `nfw_xi_t.hh` | `nfw_xi` | NFW correlation ξ(r, M) | static | `(r, M) → double` | — |

> `hmb_t.hh`/`HMB_t` (halo bias) is **retired** — bias is now read from the
> `halomodel/bias` datablock written by the Python `halo_model` module.
> The header file still exists for backwards compatibility but is not
> consumed by any active module. See Appendix A.

### 3.3 Mass-observable relations (MOR)

| Header | Type | Physics | DataBlock | `operator()` | GPU? |
|---|---|---|---|---|---|
| `mor_t.hh` | `MOR_t` | Generic Gaussian MOR P(λ\|M, z) | `cluster_abundance` (`mor_*`) | `(lt, lnM, zt)` | ✓ |
| `mor_des_t.hh` | `MOR_DES_t` | Bivariate Gaussian DES MOR | `cluster_abundance` | `(lt, lnM, zt)` | ✓ |
| `mor_des_log_t.hh` | `MOR_DES_LOG_t` | Log-space MOR (Costanzi 2019) | `cluster_abundance` | `(lt, lnM, zt)` | ✓ |
| `mor_des_2022.hh` | `MOR_DES_2022` | DES 2022 calibration | `cluster_abundance` | `(lt, lnM, zt)` | — |
| `mor_des_2022_logm.hh` | `MOR_DES_2022_LOGM` | DES 2022 with log-M parameterisation | `cluster_abundance` | `(lt, lnM, zt)` | — |
| `mor_des_triax.cuh` | `MOR_DES_TRIAX` | Triaxial-halo MOR (GPU only) | `cluster_abundance` | `(lt, lnM, zt)` | GPU |
| `mor_sdss_t.hh` | `MOR_sdss` | SDSS baseline MOR (default in `DefaultModels`) | `cluster_abundance` | `(lt, lnM, zt)` | — |
| `mor_hod_t.hh` | `MOR_HOD_t` | Shifted-Poisson HOD MOR | `cluster_mor` | `(lt, lnM, zt)` | — |
| `mor_exp_t.hh` | `MOR_EXP_t` | Exponential-tail MOR variant | `cluster_abundance` | `(lt, lnM, zt)` | — |

### 3.4 Richness / λ kernels

| Header | Type | Physics | DataBlock | `operator()` | GPU? |
|---|---|---|---|---|---|
| `lc_lt_t.hh` | `LC_LT_t` | P(λ_cen \| λ_true, z) projection + percolation (Costanzi+ 2018 eq. 15, frozen DES tables) | static `Interp2D` | `(lc, lt, zt) → double` | ✓ |
| `lc_lt_y1_t.hh` | `LC_LT_Y1_t` | P(λ_cen \| λ_true, z) Y1 form | `cluster_abundance` | `(lc, lt, zt) → double` | — |
| `lc_lt_projection_y3.hh` | `LC_LT_PROJECTION_Y3_t` | Y3 projection-only kernel | static | `(lc, lt, zt) → double` | — |
| `lo_lc_t.hh` | `LO_LC_t` | P(λ_obs \| λ_cen, R_mis) miscentering | `cluster_abundance` (`LO_LC_*`) | `(lo, lc, R_mis) → double` | ✓ |
| `int_lc_lt_des_t.hh` | `INT_LC_LT_DES_t` | Lookup-table P(λ_c \| λ_t, z) DES | static | `(lc, lt, zt) → double` | ✓ (`.cuh` `t`, `t2`, `t3`) |
| `int_lc_lt_des_t2.hh` | `INT_LC_LT_DES_t2` | Variant 2 of above | static | `(lc, lt, zt)` | ✓ |
| `int_zo_zt_t.hh` | `INT_ZO_ZT_t` | ∫P(z_obs \| z_true) dz_obs photo-z | `cluster_abundance` (`zo_zt_sigma`) | `(zomin, zomax, zt) → double` | — |
| `int_zo_zt_des_t.hh` | `INT_ZO_ZT_DES_t` | DES Gaussian photo-z (analytic) | `cluster_abundance` | `(zomin, zomax, zt) → double` | ✓ |
| `zo_zt_des_t.hh` | `ZO_ZT_DES_t` | Pointwise P(z_obs \| z_true) DES | `cluster_abundance` | `(zo, zt) → double` | — |
| `richness_kernel_t.hh` | `RICHNESS_KERNEL` | Bin-integrated K_i(λ_true, z), K_j(z) | `richness_selection` | per-bin scalars | — |
| `projection_y3_b_i_t.hh` | `PROJ_Y3_B_i_t` | B_i(λ_true, z): Y3 projection × bin integral (frozen exponential form) | static | `(lt, zt, i) → double` | — |
| `plob_ltr_emg_t.hh` | `PLOB_LTR_EMG` | EMG parameters for P(λ_obs \| λ_true, z) | `plob_ltr_params` | `(lt, zt)` | — |
| `sigma_photoz_des.hh` | `SIGMA_PHOTOZ_DES_t` | σ_phot(z_true) DES | static (z grid) | `(zt) → double` | ✓ |
| `z_kernel_data.hh` | (data) | Embedded 120-node photo-z σ(z) grid [0.10, 0.90] | static | — | — |

### 3.5 Lensing observables (NFW family)

| Header | Type | Physics | DataBlock | `operator()` | GPU? |
|---|---|---|---|---|---|
| `del_sig_t.hh` | `DEL_SIG_t` | ΔΣ Y3 1-halo | `deltasigma` | `(r, lnM, zt) → double` | — |
| `dsigma_misc.hh` | `DSIGMA_MISC` | Miscentered 1-halo NFW ΔΣ | `cluster_abundance`, `nfw_off_center` | `(r, lnM, zt)` | ✓ |
| `del_sig_sum.hh` | `DEL_SIG_SUM` | Cumulative ΔΣ | — | `(r, lnM, zt)` | — |

(2-halo Σ_prj/ΔΣ_prj/γ_t^prj live in `sigma_prj_t.hh`; see §3.7.)

### 3.6 Selection functions

| Header | Type | Physics | DataBlock | `operator()` | GPU? |
|---|---|---|---|---|---|
| `sel_function_t.hh` | `SelFunction_t` | Precomputed S_ij(lnM, z) `Interp2D` per richness × redshift bin | `sel_function` | `(i, j, lnM, zt) → double` | — |
| `pzsource_t.hh` | `PZSOURCE_t` | Binned P(z_source) | `pz_source` | `(bin, zs) → double` | — |
| `pzsource_gaussian_t.hh` | `PZSOURCE_GAUSSIAN_t` | Gaussian P(z_source) | `pz_source` | `(bin, zs) → double` | — |
| `sample_variance.hh` | `SampleVariance` | Sample-variance accumulator | — | accumulator | — |

### 3.7 Costanzi-2026 advanced operators

| Header | Type | Physics | Notes |
|---|---|---|---|
| `n_operator_sel_t.hh` | `NOperatorSelScalar<W>`, `NOperatorSelRadial<W>` | Closed-form `N_i[W]` selection-bias operator templates | Driven by all 11 `num_counts_sel` modules; `W` weight functors live in `src/modules/num_counts_sel/lensing_weights.hh`. |
| `p_operator_t.hh` | `P_operator` | Fixed-GL `P_[X]` driver: P₁, I₁, J on a common wall grid | Used by `b_sel_marg_cpu/BSelMargIntegrand.cc`. Requires `set_sample()` + `evaluate(pt)`. |
| `p_operator_cuhre_t.hh` | `P_OPERATOR_CUHRE` | Adaptive 3-D Cuhre variant of `P_operator` | Reference / convergence diagnostic. |
| `sigma_prj_t.hh` | `SIGMA_PRJ`, `ShearPrjEvaluator`, `ShearPrjCuhre`, `ShearPrjGsl` | Σ_prj + ΔΣ_prj + γ_t^prj integrand and three integrator backends | The reduced-shear `SIGMA_PRJ` form was retired 2026-05-11; the linear `ShearPrjEvaluator` is the production path. |
| `op_sel_park.hh` | `OP_SEL_PARK` | Selection operator (Park et al. form) | `.cuh` GPU variant present. |

### 3.8 Typelist

`src/models/models.hh` defines a generic bundle:

```cpp
template <typename MOR, typename LO_LC, typename LC_LT, typename ZO_ZT,
          typename ROFFSET, typename T_CEN, typename T_MIS,
          typename A_CEN, typename A_MIS, typename HMB, typename HMF,
          typename DEL_SIG, typename DV_DO_DZ, typename OMEGA_Z>
struct Models { ... };
```

`src/models/default_models.hh` instantiates the SDSS-baseline alias
`DefaultModels = Models<MOR_sdss, LO_LC_t, LC_LT_t, INT_ZO_ZT_t, ROFFSET_t,
T_CEN_t, T_MIS_t, A_CEN_t, A_MIS_t, HMB_t, HMF_t, DEL_SIG_TOM, DV_DO_DZ_t,
OMEGA_Z_SDSS>`. Most active integrands now bypass the typelist and consume
specific term types directly; the alias is preserved for older
`sigma_mort_y1` / Y1 modules. Several typelist slots (`ROFFSET`, `T_CEN`,
`T_MIS`, `A_CEN`, `A_MIS`, `HMB`, `DEL_SIG_TOM`) are retired in the
Costanzi-2026 path — see Appendix A.

### 3.9 SPT × DES cross-survey (`src/models/sptxdes/`)

| Header | Type | Physics |
|---|---|---|
| `mor_y3xspt_t.hh` | `MOR_Y3xSPT_t` | Bivariate Gaussian MOR (λ_true, ln ξ) for SPT/Y3 joint analysis |
| `mor_1d.hh` | `MOR_1D` | 1-D MOR variant |
| `omega_z_y3xspt.hh` | — | Ω(z) Y3×SPT |
| `pxizeta_t.hh` | — | P(ξ) on SPT mass-ξ parameter |
| `int_pxizeta.hh` | — | ∫P(ξ) integration |
| `ln_mez_power_law.hh` | — | Power-law ln(M)–z relation |

These headers are not currently consumed by any module in
`src/modules/CMakeLists.txt`; they are kept for the SPT cross-survey
pipeline development.

---

## 4. `src/modules` — CosmoSIS `.so` registry

Authoritative build list: [`src/modules/CMakeLists.txt`](../src/modules/CMakeLists.txt).
The order below mirrors the file's `add_subdirectory` calls.

### 4.1 Costanzi-2026 production

| Module dir | `.so` produced | Integrand template | Macro | CPU/GPU |
|---|---|---|---|---|
| `b_sel_marg_cpu` | **BSelMargIntegrand** | `y3_cluster::P_operator` | `SCALAR_EVALUATOR` | CPU fixed-GL |
| `b_sel_marg_cpu` | P1PaganiIntegrand, I1PaganiIntegrand, I2PaganiIntegrand | `NOperatorSelScalar<W>` | `SCALAR_INTEGRATION` | CPU PAGANI (benchmark only) |
| `sigma_prj_cpu` | **ShearPrjEvaluator** | `y3_cluster::ShearPrjEvaluator` | `SCALAR_EVALUATOR` | CPU fixed-GL |
| `sigma_prj_cpu` | ShearPrjCuhre | `y3_cluster::ShearPrjCuhre` | `SCALAR_EVALUATOR` | CPU adaptive (reference) |
| `sigma_prj_cpu` | ShearPrjGsl | `y3_cluster::ShearPrjGsl` | `SCALAR_EVALUATOR` | CPU GSL adaptive (diagnostic) |
| `num_counts_sel` | (see §4.2) | `NOperatorSel{Scalar,Radial}<W>` | `SCALAR_INTEGRATION` | CPU |
| `num_counts_full` | NumCountsFullScalarIntegrand | `NumCountsFullScalarIntegrand` | `SCALAR_INTEGRATION` | CPU (brute-force triple-integral reference) |

### 4.2 Selection-bias `N`-operator family (`num_counts_sel/`)

11 `.so` files, each a ≈10-line `.cc` driver typedeffing
`NOperatorSel{Scalar,Radial}<W>` with weight `W` from `lensing_weights.hh`,
all linking against `cosmosis utils models GSL CUBA fmt::fmt`.

| `.cc` file | `.so` name | Template | Computes |
|---|---|---|---|
| `NumCounts.cc` | NumCountsSel | `NOperatorSelScalar<NumCounts>` | N_i[1] selection-corrected counts |
| `MassWeighted.cc` | MassWeightedSel | `NOperatorSelScalar<MassWeighted>` | M-weighted N_i (effective mass) |
| `BiasWeighted.cc` | BiasWeightedSel | `NOperatorSelScalar<BiasWeighted>` | b-weighted N_i (effective bias) |
| `Sigma1h.cc` | Sigma1hSel | `NOperatorSelRadial<Sigma1h>` | Σ^{1h}(R) |
| `SigmaTot.cc` | SigmaTotSel | `NOperatorSelRadial<SigmaTot>` | Σ^{1h+2h}(R) |
| `DSigma1h.cc` | DSigma1hSel | `NOperatorSelRadial<DSigma1h>` | ΔΣ^{1h}(R) |
| `DSigma1hMis.cc` | DSigma1hMisSel | `NOperatorSelRadial<DSigma1hMis>` | Miscentred ΔΣ^{1h}(R) |
| `DSigmaTot.cc` | DSigmaTotSel | `NOperatorSelRadial<DSigmaTot>` | ΔΣ^{1h+2h}(R) |
| `Shear1h.cc` | Shear1hSel | `NOperatorSelRadial<Shear1h>` | γ_t^{1h}(R) |
| `Shear1hMis.cc` | Shear1hMisSel | `NOperatorSelRadial<Shear1hMis>` | Miscentred γ_t^{1h}(R) — see `pipeline_modules.tex` §11 |
| `ShearTot.cc` | ShearTotSel | `NOperatorSelRadial<ShearTot>` | γ_t^{1h+2h}(R) |

### 4.3 Validation / GPU PAGANI

| Module dir | Era / role | `.so` files | Macro |
|---|---|---|---|
| `gt_card_triax_gpu` | Triaxiality validation | `numberCounts.cu`, `avgMiscApprox{1,2}.cu`, `avgMiscFull.cu`, `avgCentFull.cu` | `CUDA_INTEGRATION` |
| `gt_mock_gpu`, `gt_mock_cpu` | Sept 2022 mock validation | numberCounts, avgCent/MiscFull, avgMiscApprox{1,2} | `CUDA` / `SCALAR_INTEGRATION` |
| `gt_card_gpu`, `gt_card_cpu` | Sept 2022 cardinal validation | same set as `gt_mock_*` | `CUDA` / `SCALAR_INTEGRATION` |
| `gt_park_sel_gpu`, `gt_park_sel_cpu` | April 2023 selection-bias σ variants | avgCentFull_sigma, avgCentFull_cum_sigma, avgMiscApprox2_sigma{,_cum} | `CUDA` / `SCALAR_INTEGRATION` |
| `sigma_buzzard_y3` | Nov 2022 Buzzard Σ/κ/γ/w_p diagnostics | avgSigmaCent/MiscBu, avgGamma{Cent,Misc,Full}Bu, avgKappa*Bu, avgWp*Bu, summedNumbersCentBu | `CUDA_INTEGRATION` |
| `sigma_park_y1` | Oct 2022 Y1 (Jim & Johnny); **currently broken**, awaits `lc_lt.cuh` debug | summedNumbers{Cent,Miscent}Y1, avgSigma{Cent,Miscent}Y1, avgWp*Y1, singleClusterPrototype | `CUDA_INTEGRATION` |

### 4.4 Y1 abundance / lensing (legacy)

| Module dir | `.so` files | Macro | Notes |
|---|---|---|---|
| `sigma_mort_y1` | sigma_cent_y1_mort, sigma_miscent_y1_mort (CPU + CUDA) | `SCALAR_INTEGRATION` / `CUDA_INTEGRATION` | Primary 2021–2022 working code; superseded by `num_counts_sel` |
| `mass_y1` | NC_cent/miscent_y1, Mass_cent/miscent_y1 | `SCALAR_INTEGRATION` | Y1 number counts and mass-weighted counts |
| `mass_mort_y1` | NC/Mass × cent/miscent (DES MOR-tality) | `SCALAR_INTEGRATION` (+ CUDA) | Y1 with DES mortality function |
| `y1_analysis` | y1_analysis, y1_analysis_mor, y1_analysis_logm_2022, y1_analysis_mor_logm_2022 | `VECTOR_INTEGRATION` | Y1 diagnostics |

### 4.5 Diagnostics

| Module dir | `.so` files | Macro | Role |
|---|---|---|---|
| `cluster_abundance_covariance` | cluster_abundance_covariance | (no integrator macro — pure module) | Covariance scaffolding |
| `snapshotsim` | Snapshotsim_ScalarIntegrand{,_NC,_max,_Sigma} | `SCALAR_INTEGRATION` | Snapshot-sim integrands |
| `buzzard_test` | buzzard_sigma_halos | `VECTOR_INTEGRATION` | Buzzard validation diagnostics |
| `sdss_analysis` | sdss_analysis | `VECTOR_INTEGRATION` | SDSS cross-checks |

### 4.6 Examples (build-system tests)

| Module dir | `.so` files | Macro |
|---|---|---|
| `ExampleScalar` | ExampleScalarIntegration_module | `SCALAR_INTEGRATION` |
| `ExampleVector` | ExampleVectorIntegration_module | `VECTOR_INTEGRATION` |
| `ExampleOneD` | ExampleOneDIntegration_module | `ONED_INTEGRATION` |

---

## 5. Module-macro primer

`src/utils/module_macros.hh` (CPU) and `src/utils/cuda_module_macros.cuh`
(GPU) define:

| Macro | When to use | Wraps |
|---|---|---|
| `DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE(klass)` | Adaptive integration of a scalar integrand class (CUBA/Cuhre under the hood). | `CosmoSISScalarIntegrationModule<klass>` |
| `DEFINE_COSMOSIS_VECTOR_INTEGRATION_MODULE(klass)` | Vector-valued integrand (multi-output). | `CosmoSISVectorIntegrationModule<klass>` |
| `DEFINE_COSMOSIS_SCALAR_EVALUATOR_MODULE(klass)` | Co-pass evaluator that does its own quadrature (e.g. fixed-GL); CosmoSIS only calls `execute()`. | `CosmoSISScalarEvaluatorModule<klass>` |
| `DEFINE_COSMOSIS_ONED_INTEGRATION_MODULE(klass)` | 1-D integration. | `OneDIntegrationModule<klass>` |
| `DEFINE_COSMOSIS_CUDA_INTEGRATION_MODULE(klass)` | GPU PAGANI integrand (in `cuda_module_macros.cuh`). | CUDA-side wrapper |

Each macro emits the C-ABI (`extern "C" setup/execute/cleanup`) that
CosmoSIS dlopens. **Exactly one** invocation per `.so`.

---

## Appendix A — Retired / not built / unused

Headers and modules that exist on disk but are not on the active path:

**Model headers (unused by any active module):**

- `hmb_t.hh` — retired 2026-05-06; bias now read from the
  `halomodel/bias` datablock written by Python `halo_model`.
- `del_sig_y1.hh`, `del_sig_tom.hh` — Y1/superseded ΔΣ forms.
- `xi_to_sigma_t.hh`, `pk_nl_t.hh`, `xinl.hh` — non-linear ξ/P(k) family;
  the active 2-halo path consumes `halomodel/{bias,xinl}` from the
  Python `halo_model` module instead.
- `roffset_t.hh` — miscentering radius prior; the Costanzi-2026 path uses
  pre-computed `nfw_dsigma_mis`/`nfw_sigma_mis` look-up tables instead.
- `t_cen_t.hh`, `t_mis_t.hh` — placeholder structs (return 1.0).
- `a_cen_t.hh`, `a_mis_t.hh` — central/satellite amplitude factors.
- `mor_exp_t.hh` — exponential MOR variant, no module instantiates it.

These slots remain in the `Models<...>` typelist so older modules
(e.g. `sigma_mort_y1`, `mass_y1`) still compile; do not delete them
without auditing those modules.

**Module directories disabled in `src/modules/CMakeLists.txt`:**

- `DESxSPTModule` — DES × SPT cross-checks (line 48).
- `sigma_kappa_y1` — Y1 convergence (line 49).
- `compton_y_sims` — Compton-y simulations (line 50).
- `model_sigmahm` — halo-model σ (line 51).

**Active dir but currently broken:**

- `sigma_park_y1` — Oct 2022 work, awaits `lc_lt.cuh` debugging.

**Orphaned subdirectories (no `CMakeLists.txt` / no `.cc`):**

- `cuda/`, `cp_camb/`, `parabola/`, `n_operator_ratios/`, `prj_lens_model/`,
  `red_shear_prj/` — stub directories under `src/modules/`.

---

## Appendix B — Python-only modules under `src/modules`

Directories that ship `.py` (+ `.ini`) but no `.cc`, so no `.so` is built:

| Dir | Role |
|---|---|
| `average_sigma_crit_inv` | ⟨Σ_crit⁻¹⟩ tabulation (consumed by `*Sel` lensing weights) |
| `deltasigma` | ΔΣ aggregation / post-processing |
| `finish` | Pipeline finalisation hooks |
| `sel_function` | S_ij(lnM, z) tabulation — see `pipeline_modules.tex` §7 |

The driving `.ini` pipelines that wire all of the above live in sibling
repos (`$PSCRATCH/github/des-cluster-nersc`,
`$PSCRATCH/github/RichnessSelection`, `$PSCRATCH/github/camb-emulator`).
