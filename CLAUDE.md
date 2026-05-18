# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Working in this tree

- **`archive/`** is a local scratch dump for files that have been retired or
  moved out of the repo but kept for reference. It's gitignored. Never write
  generated outputs, run artifacts, or new working files there — `archive/`
  is for *retiring* things, not staging new ones. If a script needs an
  output dir, write under `<ini_basename>_output/` (gitignored) or an
  explicit path the user gives you, not `archive/`.

## Project

`y3_cluster_cpp` is the DES cluster-cosmology pipeline: C++ (and CUDA) implementations
of cluster abundance, lensing, and selection-bias integrals, packaged as shared-library
**CosmoSIS modules** (`.so` files). Integration backends: CUBA/Cuhre + cubacpp (CPU)
and PAGANI (GPU) — see `externals/`, `CUBACPP_DIR`, `PAGANI_DIR`.

CosmoSIS `.ini` pipelines that drive these modules live in sibling repos
(`$PSCRATCH/github/des-cluster-nersc`, `$PSCRATCH/github/RichnessSelection`,
`$PSCRATCH/github/camb-emulator`), not in this tree. Anything still present
under `cosmosis-models/` or `y3_buzzard/` here is local/scratch and not
canonical — don't treat it as authoritative.

## Build

The project **only builds on Perlmutter GPU compute nodes with a specific non-default
toolchain**. Do not use the default module versions. Full recipe is in
`BUILD_PERLMUTTER.md`; summary:

```bash
# get a GPU node for test execution
salloc --nodes 1 --qos interactive --time 02:00:00 --constraint gpu \
       --gpus 4 --account=des_g

source ~/cosmosis_init.sh
module swap cudatoolkit/12.9 cudatoolkit/12.2   # 12.9 has nvcc bugs + missing libnvToolsExt
module swap gcc-native/13.2 gcc-native/12.3     # CUDA 12.2 rejects GCC 13+
export PATH=$(echo $PATH | tr ':' '\n' | grep -v homebrew | tr '\n' ':' | sed 's/:$//')  # homebrew pkg-config breaks Cray wrappers

export Y3_CLUSTER_CPP_DIR=/pscratch/sd/j/jesteves/y3_cluster_cpp
cd $Y3_CLUSTER_CPP_DIR && mkdir -p release-build && cd release-build
cmake -DUSE_CUDA=On -DY3GCC_TARGET_ARCH=80-real \
      -DPAGANI_DIR=${GPU_INT_DIR} -DGSL_ROOT_DIR=${CONDA_PREFIX} \
      -DCMAKE_MODULE_PATH=${CUBA_CPP_DIR}/cmake/modules \
      -DCUBACPP_DIR=${CUBA_CPP_DIR} -DCUBA_DIR=${CUBA_DIR} \
      -DCMAKE_BUILD_TYPE=Release -G Ninja ${Y3_CLUSTER_CPP_DIR}
ninja
ctest -j 10
```

- The source tree **must live on `/pscratch`**, not `/global/common/software`,
  because the latter is read-only on compute nodes.
- The bundled `externals/catch2/catch.hpp` needs the `SIGSTKSZ` patch
  (already applied in this working tree; see `BUILD_PERLMUTTER.md` §3 if it
  regresses).
- `Y3_CLUSTER_CPP_DIR` **must** point at the pscratch source root at both
  build and test time — tests use it to locate data files in `data/`.
- For CPU-only work (not recommended; most tests are GPU), omit
  `-DUSE_CUDA=On` and `Y3GCC_TARGET_ARCH`.

## Test

- Full suite: `ctest -j 10` from the build directory.
- Failing-test output: `CTEST_OUTPUT_ON_FAILURE=1 ctest`.
- Single test: run the `<name>_test` executable directly from `release-build/test/`,
  or `ctest -R <name>`. Each `*.test.cc` in `test/` compiles to its own
  executable — see `test/README.md` for how to add one.
- Non-test diagnostic executables in `test/` (e.g. `integrate_lc_lt`,
  `integrate_mor`, `n_operator_sel_t_profile`) probe individual kernels and
  accept `-h/--help`.

## CosmoSIS pipeline runs

This repo only ships the C++/CUDA modules. The `.ini` pipelines that drive
them live in sibling repos:

- `$PSCRATCH/github/RichnessSelection` — Python reference + Py↔C++ validation
  harnesses (`compare_*_py_vs_cpp.py`, `sweep_p_operator.py`).
- `$PSCRATCH/github/des-cluster-nersc` — MCMC/run-management (mock data
  vector builds, MCMC sbatch scripts, cosmology scans).
- `$PSCRATCH/github/camb-emulator` — emulator training + emulator-vs-full
  comparisons.

Run pipelines from those repos with `Y3_CLUSTER_CPP_DIR` pointing at this
tree so the built `.so` files under `release-build/src/modules/...` resolve.
Outputs land in `<ini_basename>_output/` next to the driving `.ini`.

## Architecture

### Template + macro module pattern

The core idea is that integrand *physics* lives in header-only templated
classes under `src/models/`, and a CosmoSIS module is just a `.cc` file that
instantiates the template and exports `setup/execute/cleanup` via a macro.

1. **Models** (`src/models/*.hh`, `.cuh` for device) — terms of the integrand:
   MOR (mass-observable relation), LC/LT (true/observed richness), ROFFSET
   (miscentering), HMF, DEL_SIG (ΔΣ), Ω(z), DV/dΩdz, etc. Each is a
   `struct` constructed from a `cosmosis::DataBlock`; `operator()` evaluates
   the term at the integration point.
2. **Typelist** (`src/models/models.hh`, `default_models.hh`) — packs the
   term types into a single `Models<MOR, LO_LC, LC_LT, ZO_ZT, ROFFSET,
   T_CEN, T_MIS, A_CEN, A_MIS, HMF, DEL_SIG, DV_DO_DZ, OMEGA_Z>`. Per-survey
   aliases (e.g. `DefaultModels` = SDSS set) live here. **No HMB type** —
   halo bias is now an `Interp2D` read from datablock section `haloModel/bias`
   written by `y3_buzzard/haloModelCosmosis.py`.
3. **Integrand class template** — e.g. `NOperatorSelScalar<Weight>` in
   `src/models/n_operator_sel_t.hh`. Composes the term evaluations and the
   integrator (cuhre / pagani / fixed-GL) into a thing CosmoSIS can call.
4. **Module `.cc`** (`src/modules/<name>/Foo.cc`) — `using Foo =
   SomeTemplate<Weight>;` then exactly one line:
   `DEFINE_COSMOSIS_SCALAR_INTEGRATION_MODULE(Foo);` (macro variants in
   `src/utils/module_macros.hh`: `SCALAR_INTEGRATION`, `VECTOR_INTEGRATION`,
   `SCALAR_EVALUATOR`, `ONED_INTEGRATION`, plus CUDA variants in
   `cuda_module_macros.cuh`). The macro emits the `setup/execute/cleanup`
   C-ABI that CosmoSIS dlopens.
5. **Module `CMakeLists.txt`** — one `add_library(Name MODULE Foo.cc)` +
   `target_link_libraries(... cosmosis utils models ${CUBA_LIBRARIES}
   ${GSL_LIBRARIES})`. CUDA modules add `cuda_utils cuda_models` and set
   `CUDA_ARCHITECTURES ${Y3GCC_TARGET_ARCH}`. The produced `.so` has no
   `lib` prefix (see top-level `CMakeLists.txt`).

`src/modules/CMakeLists.txt` is the registry of currently-built modules. New
modules must be added there. GPU modules are all gated on `if (USE_CUDA)`.

### Integrator conventions (Costanzi-2026 path)

- **Fixed-GL evaluator is the default**, PAGANI variants are reference-only
  benchmarks. `BSelMargIntegrand.so` co-computes P1/I1/I2 in one pass;
  `RedShearPrjEvaluator.so` co-computes Σ_prj/ΔΣ_prj/g_t^prj in one pass
  (`docs/Prj_lensing.md`).
- PAGANI variants are ~10³× slower than the fixed-GL path for these
  integrands (0.17 s vs 208 s on the b_sel wall grid) — don't treat them as
  "the production path".
- Don't vectorize scalar integrands via `scipy.integrate.quad_vec`-style
  rewrites without benchmarking: when the scalar path uses hand-placed
  panels at integrand features (e.g. `mean ± k·σ`), vectorisation destroys
  that structure (see user's global CLAUDE.md §"Performance investigation").

### Unit conventions (easy to get wrong)

These bit future-you; enforce in any new integrand:

- **chi vs R_lambda units**: `distances/d_c` returns `chi(z)` in **Mpc**.
  `R_lambda = (lam/100)^0.2` is in **cMpc/h**. Consumers must multiply
  `chi` by `h0` to put both on the same axis.
- **HMF mass shift**: `HMF_t` stores `dn/dlnM` with x-axis `ln(m_h *
  (Omega_m - Omega_nu))`. Query with `lnM + ln(Omega_m - Omega_nu)` — a
  raw `lnM` query returns the value at the wrong mass (≈ 0.6 dex shift)
  and silently produces wrong HMF values.
- **`lt` integration range**: Costanzi-2026 operators integrate lt on
  `(0, lob_centre]` *per grid point*, not a fixed `[lt_low, lt_high]`.
- **B_i and Ω(z)** are deliberately NOT in the P[X] operator — they cancel
  in every downstream b_sel ratio and are absent from the Python reference.
- **`sigma_prj` is not yet apples-to-apples with Python** — differs in NFW
  concentration (C++ c=4 vs Py c=5), `rho_s` normalisation (rho_crit vs
  rho_mean), and miscentering kernel (gamma vs delta). See
  `validations/README.md` before validating changes there.

## Python side

- `y3_buzzard/haloModelCosmosis.py` publishes `haloModel/bias` (Interp2D) and
  related halo-model outputs consumed by C++ modules. The C++ halo-bias type
  HMB_t is retired; any module that needs b(M,z) reads it from the datablock.
- Python↔C++ numeric-equivalence harnesses (`compare_*_py_vs_cpp.py`,
  `sweep_p_operator.py`) live in `$PSCRATCH/github/RichnessSelection/validations/`,
  not in this tree.

## Running on Perlmutter

Follow the NERSC rules in `~/.claude/CLAUDE.md` (shared vs debug vs preempt
QOS, `--constraint=gpu`, `--account=des_g` for GPU). Specific to this repo:

- `Y3_CLUSTER_CPP_DIR` must be in the job's env; without it tests/modules
  cannot find `data/`, `externals/`, or other sibling `.so` files.
- `.bashrc` is not sourced in batch jobs — source `cosmosis_init.sh` and
  re-do the `cudatoolkit/12.2` + `gcc-native/12.3` swaps explicitly inside
  the job script.
- GPU tests require `--gpus` on both `#SBATCH` and `srun`.
