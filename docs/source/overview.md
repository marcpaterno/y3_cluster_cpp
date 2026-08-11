# Pipeline overview and scope

This part establishes **what the repository is responsible for, what the
full CosmoSIS pipeline computes, and where the repository boundary ends**.

```{admonition} Reference configuration
:class: important
Throughout this documentation the *reference pipeline* is
`cosmosis-models/mock_mcmc_cp_camb.ini` in the **des-cluster-nersc**
repository (commit `5055375`, 2026-08-09): the production mock-MCMC
configuration, wiring the miscentering-aware one-halo branch
(`Shear1hMisSel`) and the fixed-GL projection evaluator
(`ShearPrjEvaluator`, section `shear_prj`). The copy of the same file
inside `y3_cluster_cpp/cosmosis-models/` is its single-sample smoke
variant (test sampler, centered-only `Shear1hSel`) and is documented as
such.
```

## Purpose of `y3_cluster_cpp`

`y3_cluster_cpp` provides the high-performance numerical implementation of
the cluster-population and cluster-lensing calculations used in the DES
cluster-cosmology pipeline.

The repository is organised as a collection of CosmoSIS modules. These
modules read cosmological, halo-model, richness–mass relation,
survey-selection, and lensing-calibration quantities from the CosmoSIS
DataBlock; evaluate the required population integrals; and write predicted
cluster observables back to the DataBlock.

The main calculations implemented by the current Costanzi-2026 path are:

1. expected cluster number counts in observed-richness and observed-redshift
   bins;
2. the population-averaged one-halo tangential shear;
3. the lensing contribution from line-of-sight structures associated with
   optical projection selection;
4. intermediate quantities required to predict the scale-dependent optical
   selection bias;
5. diagnostic quantities such as effective mass, effective halo bias,
   $\Sigma_{\rm prj}$, and $\Delta\Sigma_{\rm prj}$.

The repository is therefore best understood as the **cluster-observable
prediction engine**. The complete cosmological analysis also depends on
CosmoSIS Standard Library modules, local Python modules, data and
calibration products, and external run-management configurations. The
documentation distinguishes four things that are easily conflated:

- the **scientific pipeline** (the sequence of calculations);
- the **modules implemented by this repository**;
- the **external CosmoSIS and Python modules**;
- the **driver configuration used to assemble them**.

## Scientific outputs

The central outputs are the expected cluster abundance and stacked
tangential shear in bins of observed richness $\lambda_{\rm ob}$ and
observed redshift $z_{\rm ob}$. The reference binning is 4 richness bins,
$\lambda_{\rm ob} \in \{[20,30),\,[30,45),\,[45,60),\,[60,200)\}$, times 3
redshift bins, $z_{\rm ob} \in \{[0.20,0.35),\,[0.35,0.50),\,[0.50,0.65)\}$
— 12 bins in total, with 10 radial points per bin for the lensing profile
(a 12-point number-count vector and a 120-point shear vector).

### Number counts

For observed-richness bin $i$ and observed-redshift bin $j$, the
number-count calculation has the general form

$$N_{ij} = \int dz_{\rm true} \int d\ln M\;
\Omega(z_{\rm true})\,
\frac{dV}{d\Omega\,dz_{\rm true}}\,
\frac{dn}{d\ln M}(M, z_{\rm true})\,
S_{ij}(M, z_{\rm true}),$$

where

- $\Omega(z)$ is the effective survey area;
- $dV/d\Omega\,dz$ is the comoving volume element;
- $dn/d\ln M$ is the halo mass function;
- $S_{ij}(M, z)$ is the joint richness and redshift selection probability.

The selection probability connects the halo population to the observed
cluster catalogue. Its richness component is ultimately based on

$$P(\lambda_{\rm ob}\mid M, z) = \int d\lambda_{\rm true}\,
P(\lambda_{\rm ob}\mid\lambda_{\rm true}, z)\,
P(\lambda_{\rm true}\mid M, z).$$

This separates the intrinsic richness–mass relation from observational
scatter, projection effects, and related cluster-finder effects. The
factorisation follows the Costanzi cluster-selection framework used by the
DES analyses.

### Population-weighted one-halo lensing

The one-halo module does not directly return a single representative
cluster profile. It evaluates a population integral,

$$N_{ij}\!\left[\gamma_t^{1h}\right](R) = \int dz \int d\ln M\;
\Omega(z)\,\frac{dV}{d\Omega\,dz}\,\frac{dn}{d\ln M}\,
S_{ij}(M, z)\,\gamma_t^{1h}(R; M, z).$$

The stacked, per-cluster mean profile is obtained by dividing by the
corresponding number-count prediction:

$$\left\langle\gamma_t^{1h}(R)\right\rangle_{ij}
= \frac{N_{ij}\!\left[\gamma_t^{1h}\right](R)}{N_{ij}[1]}.$$

In the reference (miscentering-aware) implementation, `Shear1hMisSel`, the
one-halo weight is the centered-plus-miscentered mixture

$$\Delta\Sigma_{1h}^{\rm full}
= (1 - f_{\rm mis})\,\Delta\Sigma_{\rm cen}
+ f_{\rm mis}\,\Delta\Sigma_{\rm mis},$$

with $(f_{\rm mis}, \tau_{\rm mis})$ read from the DataBlock
`miscentering` section when present. `Shear1hSel` is the centered-only
variant ($f_{\rm mis}=0$), used by the smoke pipeline. The corresponding
shear is obtained by multiplying by the source-weighted critical
surface-density inverse.

### Projection-selected lensing contribution

The projection branch predicts the contribution to the stacked density and
shear profiles from foreground and background haloes lying along the
cluster line of sight.

This contribution is physically distinct from the ordinary one-halo
population integral. It depends on:

- the halo mass function;
- halo bias;
- the nonlinear matter correlation function;
- the redshift projection kernel;
- the angular overlap between the target cluster and projected systems;
- the scale-dependent optical selection bias.

The module produces

$$\Sigma_{\rm prj}(R), \qquad \Delta\Sigma_{\rm prj}(R), \qquad
\gamma_t^{\rm prj}(R).$$

The projection result is already evaluated per observed-richness and
observed-redshift bin, and therefore is **not** divided by `NumCountsSel`.

### Final shear prediction

For the current tangential-shear formulation, the theoretical stacked
profile is assembled as

$$\boxed{\;\gamma_t^{\rm theory}(R \mid i, j)
= \frac{N_{ij}[\gamma_t^{1h,\rm full}](R)}{N_{ij}[1]}
+ \gamma_t^{\rm prj}(R \mid i, j)\;}$$

(with $\gamma_t^{1h,\rm full} \to \gamma_t^{1h}$ in the centered-only
variant). The two terms can be added directly because both are predictions
for tangential shear $\gamma_t$, **not** reduced shear $g_t$: the retired
reduced-shear formulation $g_t = \gamma_t/(1-\kappa)$ has a denominator
that prevents this additive decomposition (see
{doc}`systematics/index`).

```{note}
The reference mock configuration runs `average_sigma_crit_inv` with
`unity = T`, i.e. $\langle\Sigma_{\rm crit}^{-1}\rangle \equiv 1$, so the
quantity actually emitted and fitted is $\Delta\Sigma(R)$. This keeps the
mock data vector and the theory prediction the same observable. With a
physical $\langle\Sigma_{\rm crit}^{-1}\rangle(z_l)$ the identical
composition yields $\gamma_t(R)$.
```

## Pipeline boundary

The complete calculation contains three categories of modules.

### Upstream cosmology and halo-population modules

These generate the cosmological quantities consumed by the cluster
calculations:

- cosmological parameter consistency;
- linear growth;
- linear and nonlinear matter power spectra;
- cosmological distances;
- the halo mass function;
- halo bias;
- the nonlinear matter correlation function.

Some are supplied by the CosmoSIS Standard Library, others are local Python
or emulator modules.

### Cluster-specific modules

These implement the DES cluster model:

- richness and redshift selection;
- cluster number counts;
- one-halo lensing with miscentering;
- optical selection-bias operators;
- the projected density and shear contribution;
- derived population-weighted quantities.

The computationally expensive population and projection integrals are
primarily implemented in C++, with CUDA or alternative integration backends
retained for selected modules, benchmarks, or historical implementations.

### Likelihood and run management

The likelihood module assembles the final observable vector and compares it
with a data vector and covariance matrix. MCMC configuration, job
submission, chain management, and Python-to-C++ validation workflows live
outside the core compiled-module layer, in sibling repositories
(`des-cluster-nersc`, `RichnessSelection`, `camb-emulator`).

No single `.ini` file inside `y3_cluster_cpp` uniquely defines the
production DES analysis; the reference configuration named at the top of
this page, with its commit, is the traced example.

## High-level data flow

```text
Cosmological parameters
HOD / richness–mass parameters
Survey and nuisance parameters
                |
                v
   Cosmological consistency and growth
                |
                v
      Matter power spectrum P(k,z)
                |
                v
   Halo mass function and halo model
      n(M,z), b(M,z), xi_NL(r,z)
                |
        +-------+-------------------+
        |                           |
        v                           v
Richness/redshift selection   Lensing geometry
      S_ij(M,z)               <Sigma_crit^-1>(z)
        |                           |
        +------------+--------------+
                     |
          +----------+-----------+
          |                      |
          v                      v
   Number-count and       Optical selection-bias
   one-halo operators      P1, I1, J -> b_sel
          |                      |
          |                      v
          |             Projection density profile
          |             Sigma_prj, DeltaSigma_prj
          |                      |
          v                      v
      N_ij and             gamma_t^prj(R)
 N_ij[gamma_t^1h](R)              |
          |                       |
          +-----------+-----------+
                      |
                      v
             Final theory vector
       N_ij and gamma_t^theory(R)
                      |
                      v
                  Likelihood
```

The reference CosmoSIS configuration follows this same sequence: cosmology,
power spectrum, mass function, halo model, lensing geometry, selection,
counts and one-halo shear, selection-bias operators, projection shear, and
likelihood.

## Module map

The table below describes **the reference Costanzi-2026 path**, not every
directory registered in `src/modules/CMakeLists.txt`. The module registry
still contains historical Y1 implementations, mock-validation modules,
CPU/GPU experiments, diagnostic backends, and examples — see the
{doc}`historical appendix <modules/historical>`.

| Pipeline layer | Modules (ini section) | Main products |
|---|---|---|
| Cosmological consistency | `consistency`, `GrowthFactor` | consistent parameters, $D(z)$ |
| Power spectrum | `cp_camb` (CosmoPower emulator; CAMB is the full-fidelity path) | $P_{\rm lin}(k,z)$ and related products |
| Halo population | `MfTinker`, `halo_model` | $dn/d\ln M$, $b(M,z)$, $\xi_{\rm NL}(r,z)$, NFW profiles |
| Projection calibration | `y3_buzzard/prj_params.py` (imported directly by `bsel.py`; retired as a standalone module) | parameters of $P(\lambda_{\rm ob}\mid\lambda_{\rm true},z)$ |
| Lensing geometry | `average_sigma_crit_inv` | $\langle\Sigma_{\rm crit}^{-1}\rangle(z_l)$ |
| Selection | `sel_function` | $S_{ij}(\ln M, z)$ |
| Cluster counts | `NumCountsSel` | $N_{ij}[1]$ |
| One-halo lensing | `Shear1hMisSel` (reference) / `Shear1hSel` (centered-only) | $N_{ij}[\gamma_t^{1h}](R)$ |
| Selection-bias operators | `b_sel_marg`, `bsel` | $P_1$, $I_1$, $J$, and scale-dependent $b_{\rm sel}$ |
| Projection lensing | `shear_prj` (`ShearPrjEvaluator`) | $\Sigma_{\rm prj}$, $\Delta\Sigma_{\rm prj}$, $\gamma_t^{\rm prj}$ |
| Likelihood | `y3_buzzard/likelihood_cp.py` | final log-likelihood |

## Repository architecture

The active C++ code follows a separation between physical model terms and
executable CosmoSIS modules. This model–module separation is a central
design principle; it is introduced here and used throughout the
module-by-module documentation.

**`src/models/`** contains reusable physical and numerical components,
generally implemented as header-only templates:

- the richness–mass relation;
- observed-versus-true richness kernels;
- photometric-redshift kernels;
- the halo mass function;
- volume and survey-area factors;
- NFW and miscentered lensing profiles;
- generic population-integral operators.

Each model term evaluates one component of an integrand using information
read from the CosmoSIS DataBlock.

**`src/modules/`** contains the CosmoSIS-facing module implementations. A
module generally:

1. selects a model or weight;
2. instantiates the appropriate integration template;
3. selects an integration backend;
4. exposes the CosmoSIS `setup`, `execute`, and `cleanup` interface;
5. writes its outputs to a named DataBlock section.

**`src/utils/`** contains shared infrastructure: one- and two-dimensional
interpolation, DataBlock readers, numerical integration utilities, CosmoSIS
module macros, and CPU/CUDA support utilities.

## The DataBlock as the pipeline interface

The CosmoSIS DataBlock is the contract connecting all pipeline stages, and
the only supported communication mechanism between modules. Each module is
documented in terms of

$$\boxed{\text{DataBlock inputs} \;\rightarrow\; \text{calculation}
\;\rightarrow\; \text{DataBlock outputs}}$$

For example:

```text
sel_function
    reads:  cosmological parameters; richness–mass parameters;
            projection-kernel parameters; richness and redshift bins
    writes: selection/S_stack = S_ij(lnM, z)

NumCountsSel
    reads:  selection/S_stack; mass_function/*; distances/*; survey area
    writes: numcountssel/vals

Shear1hMisSel  (or Shear1hSel)
    reads:  selection/S_stack; haloModel/dSigma_nfw;
            average_sigma_crit_inv/*; miscentering parameters (optional)
    writes: shear1hmissel/vals  (or shear1hsel/vals)

shear_prj
    reads:  mass function; halo bias; nonlinear correlation function;
            optical selection-bias products; projection kernel;
            lensing geometry
    writes: sigma_prj/*, dsigma_prj/*, shear_prj/*
```

The detailed key names and units belong to the
{doc}`module reference <modules/index>`.

## Terminology

To avoid ambiguity, the following terms are fixed here and used throughout.

| Term | Symbol | Definition |
|---|---|---|
| True richness | $\lambda_{\rm true}$ | The intrinsic cluster richness before observational scatter and projection contamination. |
| Observed richness | $\lambda_{\rm ob}$ | The richness assigned by the cluster finder and used to select the observed sample. |
| One-halo lensing | — | The projected density or shear generated by the main halo associated with the selected cluster. |
| Standard two-halo term | — | The conventional correlated-matter contribution around a halo, weighted by its halo bias. |
| Projection-selected lensing contribution | — | The contribution from line-of-sight haloes whose presence also affects the optical selection. This is the Costanzi-2026 `shear_prj` branch; it should **not** be referred to as the ordinary two-halo term. |
| Tangential shear | $\gamma_t$ | The current primary lensing observable, $\gamma_t = \Delta\Sigma\,\langle\Sigma_{\rm crit}^{-1}\rangle$. |
| Reduced shear | $g_t$ | The nonlinear observable $g_t = \gamma_t/(1-\kappa)$. **Not used by the current reference path** (retired; see {doc}`systematics/index`). |
| Lensing boost factor | — | A correction for contamination of the source sample by cluster-associated galaxies. Distinct from both $\Sigma_{\rm crit}^{-1}$ (lensing geometry) and the projection-selection bias. No boost-factor module exists in the reference pipeline. |

## Status decisions

The following decisions were fixed on 2026-08-10 and are reflected
throughout this documentation:

1. **Canonical configuration**: `des-cluster-nersc
   cosmosis-models/mock_mcmc_cp_camb.ini` @ `5055375`; the in-repo copy is
   its smoke variant.
2. **Reference one-halo branch**: `Shear1hMisSel` (miscentering-aware);
   `Shear1hSel` is the centered-only variant.
3. **Naming**: the projection stage is called `shear_prj` (the DataBlock
   output section); `red_shear_prj` is a legacy name for the same branch.
4. **Data-vector length**: 12 number counts + 12 bins × 10 radii = 120
   shear points (`likelihood_cp.py` `_SHEAR_N = 120`; older comments
   quoting 180/15-radii are stale).
5. **Reduced shear**: retired/historical; tangential shear (or ΔΣ under
   the unity convention) is the primary observable.
