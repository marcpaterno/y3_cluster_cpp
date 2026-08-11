# Validation and testing

## Existing baseline

What already exists, mapped to what it covers:

- **Unit tests** — 44 Catch2 `*.test.cc` executables plus 7 CUDA tests
  (`cuda_interp_1d/2d`, `nfw_sigma_mis`, `nfw_dsigma_mis`,
  `int_lc_lt_des_t/f`, `timing_sentry`), run via `ctest` from the build
  directory. Includes the `*_probability_unity` normalisation checks
  (lc_lt, lo_lc, mor, roffset).
- **Diagnostic executables** — `integrate_delta_sigma`,
  `integrate_lc_lt`, `integrate_lo_lc`, `integrate_mor`,
  `integrate_roffset` probe individual kernels (`-h/--help`).
- **Python ↔ C++ harnesses** — `compare_*_py_vs_cpp.py` and
  `sweep_p_operator.py` in `RichnessSelection/validations/`.
- **Emulator validation** — `docs/emulator_validation.tex`: `cp_camb`
  vs CAMB on number counts at the Planck fiducial.

## The closure test (headline end-to-end check)

The reference configuration doubles as the regression test:

1. run the test sampler at the fiducial parameters — the Gaussian
   likelihood against the mock data vector must give
   $\log L \approx 0$ (tolerance $|\log L| < 10^{-2}$; the
   Ω(z)-convention experiment showed how a convention mismatch breaks
   this to $-151.7$);
2. run the production sampler — posteriors must recover the fiducial
   HOD parameters within their 68% intervals.

**Executed 2026-08-10 on a Perlmutter login node: the closure test
passes with $\log L = 0.0$.**

## Full Pipeline Recipe

The reference regression recipe, as executed (all CPU, single thread,
$\sim 1.5$ s per pipeline evaluation; configurations in
`des-cluster-nersc/cosmosis-models/`):

1. **Environment**: `source ~/cosmosis_init.sh`; export
   `DES_CLUSTER_NERSC_DIR` and put `$Y3_CLUSTER_CPP_DIR` on
   `PYTHONPATH` (needed for the `y3_buzzard` package imports in
   `bsel.py`/`likelihood_cp.py`).
2. **Generate the mock data vector**:
   `cosmosis generate_mock_dv_shear_prj.ini` — the full pipeline at the
   fiducial point (test sampler) plus the packaging module; writes
   `$Y3_CLUSTER_CPP_DIR/data/mock/mock_dv_cp_camb.npz` (12 counts +
   120 ΔΣ points; Poisson NC inverse covariance; Buzzard jackknife
   shear covariance rebinned 20 → 10 radii).
3. **Closure**: `cosmosis mock_closure_test.ini` — an independent
   pipeline evaluation plus `likelihood_cp.py` against the stored
   vector. Expected: `Likelihood = 0.0`.

**Fiducial point** (`mock_mcmc_widePlanck_values.ini`): $h_0 = 0.6766$,
$\Omega_m = 0.311049$, $\Omega_b = 0.048975$, $n_s = 0.9665$,
$\sigma_8 = 0.8238$, $m_\nu = 0$; HOD
$\log_{10} M_{\rm min} = 11.4$, $\log_{10}(M_1/M_{\rm min}) = 1.3$,
$\alpha = 0.86$, $\epsilon = 0$, $\sigma_\lambda = 0.18$; `hmf_s` = 0,
`hmf_q` = 1; $f_{\rm mis} = 0.22$, $\tau_{\rm mis} = 0.17$;
`unity = T` (observable is ΔΣ).

**Expected number counts** $N_i[1]$ (richness-fast, z-block-slow;
regression values from the executed run):

| | $\lambda\in[20,30)$ | $[30,45)$ | $[45,60)$ | $[60,200)$ |
|---|---:|---:|---:|---:|
| $z\in[0.20,0.35)$ | 1380.58 | 513.75 | 141.43 | 91.72 |
| $z\in[0.35,0.50)$ | 2300.09 | 784.91 | 196.28 | 110.29 |
| $z\in[0.50,0.65)$ | 2473.19 | 778.78 | 177.74 | 87.88 |

**Expected lensing prediction**: the 120-point ΔΣ vector stored as
`data_Shear` in the generated `.npz` (range 7.49–148.0
$h\,M_\odot/\mathrm{pc}^2$ at fiducial); its two-composition breakdown
is shown in the {doc}`science chapter figure <../science/index>`.
Regenerating the mock and rerunning the closure after any code change
is the end-to-end regression check.

## Test policy for new numerical routines

Each important numerical routine should ideally have:

- a simple **analytic** test where possible;
- a **high-accuracy reference** numerical calculation;
- a **regression** test against the production implementation.

Which of the three applies is decided per routine, case by case, before
implementation.

## CPU–GPU validation

```{todo}
For every GPU implementation: CPU ↔ GPU ↔ high-precision reference.
Record fractional differences, tolerances, runtime, and hardware.
Where a GPU variant is broken or retired (e.g. `sigma_park_y1`), record
that status instead of implying coverage.
```

## Known open discrepancies

Documented, not hidden:

- **`sigma_prj` C++ vs Python** is not yet apples-to-apples: NFW
  concentration $c=4$ (C++) vs $c=5$ (Python reference); $\rho_s$
  normalisation (rho_crit vs rho_mean); miscentering kernel (gamma vs
  delta).
- **Python HMF z-evolution**: the Python reference must keep the
  Tinker-08 $(1+z)$ scaling enabled to match C++ (closes an ~8% residual
  to <1%).
- **Stale comments**: several ini/docstring comments still quote a
  180-point (15-radii) shear vector; the production vector is
  12 × 10 = 120 (`likelihood_cp.py`).
