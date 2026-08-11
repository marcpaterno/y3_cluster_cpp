# Systematic models

## Miscentering

Miscentering enters the model in two distinct places:

- miscentering of the **lensing profile** — the reference one-halo branch
  `Shear1hMisSel` uses the DES-Y3 redMaPPer offset model (Kelly et al.
  2023 gamma kernel), mixing centered and miscentered profiles with
  $(f_{\rm mis}, \tau_{\rm mis})$ read from the `miscentering` DataBlock
  section (defaults $f_{\rm mis}=0.22$, $\tau_{\rm mis}=0.17$). The full
  model — offset kernels, azimuthal convolution, and the precomputed
  lookup tables — is documented in
  {doc}`../science/index` §One-halo lensing and miscentering.
- miscentering of the **abundance selection**.

```{todo}
Document the abundance-selection side of miscentering, and state
explicitly which effects are forward modeled and which enter through
precomputed corrections.
```

## Optical selection bias

The mean lensing surface density around optically selected clusters of
observed richness $\lambda^{\rm ob}$ at observed redshift $z^{\rm ob}$
splits into a one-halo term (the cluster's own miscentered NFW profile) and
a two-halo term from correlated neighbouring halos along the line of sight,

$$
\langle\Sigma(R\mid\lambda^{\rm ob},z^{\rm ob})\rangle \;=\;
\langle\Sigma^{\rm 1h}(R\mid\lambda^{\rm ob},z^{\rm ob})\rangle \;+\;
\langle\Sigma^{\rm prj}(R\mid\lambda^{\rm ob},z^{\rm ob})\rangle .
$$

Because redMaPPer counts richness through red-sequence galaxies — biased
tracers of large-scale structure — a line of sight with more red galaxies
gets assigned a larger $\lambda^{\rm ob}$ at fixed halo mass. The two-halo
term at fixed $\lambda^{\rm ob}$ therefore does not carry the bare halo
bias $b(M,z)$: it carries a *selection-affected* cluster bias
$b_{\rm sel}(\lambda^{\rm ob},z^{\rm ob},\theta)$ that interpolates between
a small-scale limit dominated by close red-galaxy neighbours,
$b_{\rm sel}^{\rm small}$, and a large-scale limit set by correlated LSS,
$b_{\rm sel}^{\rm large}$.

### Halo-bias aggregate

The anchor of the pipeline is a single 1-D $\ln M$ integral at fixed
$z^{\rm ob}$, weighting the Tinker-2010 halo bias by the halo-mass
distribution at fixed observed richness:

$$
b_{\rm halo}(\lambda^{\rm ob},z^{\rm ob})
 = \frac{\int d\ln M\;M\,n(M,z^{\rm ob})\,b(M,z^{\rm ob})\,P(\lambda^{\rm ob}\mid M,z^{\rm ob})}
        {\int d\ln M\;M\,n(M,z^{\rm ob})\,P(\lambda^{\rm ob}\mid M,z^{\rm ob})} ,
$$

with $n(M,z)$ the Tinker-2008 halo mass function, $b(M,z)$ the Tinker-2010
peak-height bias, $P(\lambda\mid M,z)$ the mass–observable relation pdf,
and $\ln M$ running over $[10^{13},10^{15.5}]\,h^{-1}M_\odot$.

### Sigmoid ansatz

The $\theta$-shape of the selection bias is fixed analytically (Costanzi
et al. 2026); only the two plateaus need computing:

$$
b_{\rm sel}(\lambda^{\rm ob},\lambda^{\rm tr},z^{\rm ob},\theta)
= b_{\rm sel}^{\rm small}\,[\,1-\sigma(\theta)\,] + b_{\rm sel}^{\rm large}\,\sigma(\theta),
\qquad
\sigma(\theta) = \frac{1}{1+e^{-k(\theta-\theta_0)}},
$$

with $k=2.5/\theta_\lambda$, $\theta_0=\theta_\lambda/2$, and
$\theta_\lambda \equiv R_\lambda(\lambda^{\rm ob})/D_A(z^{\rm ob})$ the
angle subtended by the target's redMaPPer radius
$R_\lambda(\lambda^{\rm ob})=(\lambda^{\rm ob}/100)^{0.2}\,h^{-1}\mathrm{Mpc}$.

### The $\mathcal{P}[X]$ operators

The projection kernel is the per-halo weight with which a projected halo
at $(\lambda,z,\theta)$ contaminates the target's observed richness,

$$
\rho^{\rm prj}(\lambda,z,\theta\mid\lambda^{\rm ob},z^{\rm ob})
= w_z(z,z^{\rm ob})\,f_A(\theta,\lambda^{\rm ob},z^{\rm ob},\lambda,z)\,\lambda,
\qquad 0<\lambda<\lambda^{\rm ob},
$$

and $\rho^{\rm prj}\equiv 0$ for $\lambda\ge\lambda^{\rm ob}$ — the
redMaPPer ranking cut: only lower-richness projected halos contaminate the
target, so every $\lambda$-integral runs to upper limit $\lambda^{\rm ob}$.
Here $w_z$ is the parabolic photo-$z$ kernel and $f_A\in[0,1]$ is the
closed-form overlap fraction of the two redMaPPer disks,
$f_A = A_{\rm ov}(\theta,\theta_{\lambda,{\rm ob}},\theta_\lambda(\lambda,z))
/[\pi\,\theta_\lambda^2(\lambda,z)]$ — normalised by the *projector's* own
disk area, with $A_{\rm ov}$ the circle–circle intersection area at centre
separation $\theta$. $f_A$ vanishes beyond
$\theta_{\lambda,{\rm ob}}+\theta_\lambda$, so the $\theta$-support is
finite.

The $\mathcal{P}[X]$ operator is the line-of-sight-averaged integral of an
arbitrary $X(M,z,\theta)$ weighted by the halo distribution and
$\rho^{\rm prj}$:

$$
\mathcal{P}[X](\lambda^{\rm ob},z^{\rm ob})
\equiv
\int\!dz\,\frac{dV}{dz\,d\Omega}(z)\!
\int\!dM\,n(M,z)\!
\int\!d\lambda\,P(\lambda\mid M,z)\;
2\pi\!\int\!d\theta\,\sin\theta\;
\rho^{\rm prj}(\lambda,z,\theta\mid\lambda^{\rm ob},z^{\rm ob})\,X(M,z,\theta\mid z^{\rm ob}).
$$

Three specialisations drive the pipeline:

$$
\mathcal{P}[1] = \langle\Delta^{\rm prj}_{\rm bkg}\rangle(\lambda^{\rm ob},z^{\rm ob}),
\qquad
I_2 \equiv \mathcal{P}\!\left[b\,\xi_{\rm NL}\right],
\qquad
I_1 \equiv \mathcal{P}\!\left[b\,\xi_{\rm NL}\,\sigma(\theta)\right],
$$

where $\xi_{\rm NL}(z,z^{\rm ob},\theta)$ is the nonlinear matter
correlation function at the 3-D comoving separation between target and
projected halo; by linearity
$\mathcal{P}[b\,\xi_{\rm NL}(1-\sigma)] = I_2 - I_1$. Halo exclusion sets
$\xi_{\rm NL}\equiv 0$ whenever $\Delta\chi<R_{\rm excl}$, with
$R_{\rm excl}\equiv R_\lambda(\lambda^{\rm ob})(1+z^{\rm ob})$.

The plateaus then close in two steps. The random-aperture contamination is
$\Delta^{\rm prj}_{\rm RND} = \mathcal{P}[1] + b_{\rm halo}\,I_2$; the
large-$\theta$ plateau follows from the Buzzard-calibrated relation

$$
b_{\rm sel}^{\rm large}(\lambda^{\rm ob},\lambda^{\rm tr},z^{\rm ob}) = b_{\rm halo}\,
  \bigl[\,1 + 0.13\,\delta^{\rm prj}\,\bigr],
\qquad
\delta^{\rm prj} \equiv \frac{\lambda^{\rm ob}-\lambda^{\rm tr}}{\Delta^{\rm prj}_{\rm RND}} - 1,
$$

and the small-$\theta$ plateau closes the remaining linear equation:

$$
b_{\rm sel}^{\rm small}(\lambda^{\rm ob},\lambda^{\rm tr},z^{\rm ob}) =
\frac{(\lambda^{\rm ob}-\lambda^{\rm tr}) - \mathcal{P}[1] - b_{\rm sel}^{\rm large}\,I_1}{I_2 - I_1}.
$$

### Marginalisation over $\lambda^{\rm tr}$

The final input to the lensing profile is the true-richness marginal,

$$
b_{\rm sel}(\lambda^{\rm ob},z^{\rm ob},\theta) =
\int\!d\lambda^{\rm tr}\;b_{\rm sel}(\lambda^{\rm ob},\lambda^{\rm tr},z^{\rm ob},\theta)\,
P(\lambda^{\rm tr}\mid\lambda^{\rm ob},z^{\rm ob}),
$$

with $P(\lambda^{\rm tr}\mid\lambda^{\rm ob},z^{\rm ob})\propto
P(\lambda^{\rm ob}\mid\lambda^{\rm tr},z^{\rm ob})\,\pi(\lambda^{\rm tr},z^{\rm ob})$
and $\pi(\lambda^{\rm tr},z^{\rm ob}) = \int dM\,n(M,z)\,
P(\lambda^{\rm tr}\mid M,z)$ the true-richness prior. Because the sigmoid
carries no $\lambda^{\rm tr}$ dependence, the marginalisation commutes
with the sigmoid assembly: defining the $\lambda^{\rm tr}$-averaged
plateaus $\bar b_{\rm small}$, $\bar b_{\rm large}$ (same average applied
to each plateau),

$$
b_{\rm sel}(\lambda^{\rm ob},z^{\rm ob},\theta)
= \bar{b}_{\rm small}(\lambda^{\rm ob},z^{\rm ob})\,[\,1-\sigma(\theta)\,]
+ \bar{b}_{\rm large}(\lambda^{\rm ob},z^{\rm ob})\,\sigma(\theta),
$$

an exact algebraic rearrangement that moves the
$\lambda^{\rm tr}$-quadrature out of the $\theta$-loop, so downstream
consumers pay only the analytic sigmoid per call.

### Numerical splits

The three forward integrals $\mathcal{P}[1]$, $I_1$, $I_2$ are 4-D
integrals over $(z,\theta,M,\lambda)$. The inner $(M,\lambda)$ axes
converge with plain Gauss–Legendre ($N_M=24$, $N_\lambda=60$); the outer
axes carry two sharp, physics-driven features:

- **$\theta$-axis: split at the exclusion scale.** At each $z$, halo
  exclusion kills $\theta$ below $\theta_{\rm excl}(z)$, defined by
  $\cos\theta_{\rm excl}(z) =
  [\chi(z)^2+\chi(z^{\rm ob})^2-R_{\rm excl}^2]/[2\,\chi(z)\,\chi(z^{\rm ob})]$
  (clipped to $[-1,1]$). A fixed GL grid on $(0,2\theta_\lambda)$ with a
  mask integrates across this discontinuity and stalls at a $\sim 0.3\%$
  error floor. The fix is a per-$z$ change of interval: integrate on
  $[\theta_{\rm excl}(z),\,2\theta_\lambda]$ with no mask, so the step
  sits at an endpoint where GL nodes cluster. $N_\theta=10$ then reaches
  sub-$0.01\%$ — a $\sim 12\times$ node reduction and two decades better
  precision than the mask recipe.
- **$z$-axis: split by physics.** Each operator factorises into an outer
  $z$-integral over
  $f_X(z) = (dV/dz\,d\Omega)\,w_z(z,z^{\rm ob})\,\mathcal{I}^X_{\rm inner}(z)$.
  For $I_1,I_2$ the exclusion cut carves three regions around
  $z^{\rm ob}$, with
  $\delta z_{\rm excl}\equiv R_{\rm excl}/(c/H(z^{\rm ob}))\approx
  4.7\times10^{-4}$ at the reference point: a flat *ring plateau* at
  $|z-z^{\rm ob}|<\delta z_{\rm excl}$ (only $\Delta\chi>R_{\rm excl}$
  structure contributes), twin *exclusion peaks* at
  $z=z^{\rm ob}\pm\delta z_{\rm excl}$ (where $\xi_{\rm NL}$ is evaluated
  at $\Delta\chi\to R_{\rm excl}$, its largest resolved value), and an
  *outer decay* $\propto\Delta\chi_\parallel^{-1.7}$ modulated by $w_z$
  out to $|z-z^{\rm ob}|\sim\sigma_z(z^{\rm ob})\approx 0.09$. The
  integral is written as
  $I = I_{\rm ring} + I^{\rm fg}_{\rm outer} + I^{\rm bg}_{\rm outer}$:
  GL in $z$ on the ring band, GL in $u=\ln|\Delta\chi_\parallel|$ on the
  two outer halves so the $\xi_{\rm NL}\propto e^{-1.7u}$ peak lands at
  the inner endpoint. Sharp features sit at region boundaries, and at
  $N_z=80$ all three operators reach sub-$0.01\%$ against
  `scipy.integrate.quad` with matched inner grids.

*Source: `richness_selection.tex` §Scope and notation, §Core equations,
§The θ-axis: split at exclusion, §The z-axis: split by physics, §Glossary
of quantities.*

## Lensing geometry and boost factors

Two concepts that must not be conflated:

**Lensing geometry** — the `average_sigma_crit_inv` module computes the
source-distribution-weighted lensing efficiency
$\langle\Sigma_{\rm crit}^{-1}\rangle(z_l)$. This converts excess surface
density into tangential shear,
$\gamma_t = \Delta\Sigma\,\langle\Sigma_{\rm crit}^{-1}\rangle$. Under the
`unity = T` convention of the reference mock configuration it is set to 1
and the pipeline observable is $\Delta\Sigma$ itself.

**Lensing boost factor** — a correction for contamination of the *source*
sample by cluster-associated galaxies (radially dependent, largest at
small $R$). This is a different physical effect from both
$\Sigma_{\rm crit}^{-1}$ and the projection-selection bias.
**No boost-factor module exists in the reference pipeline**; this section
documents the concept and where such a correction would enter, so the
three notions stay distinct.

```{todo}
Document the boost-factor model (physical origin, radial dependence,
relation to $\Sigma_{\rm crit}^{-1}$, and whether a correction would
apply to shear, $\Delta\Sigma$, or another intermediate quantity) if and
when it is added to the analysis.
```

## Tangential vs reduced shear

The primary lensing observable of the current path is the **tangential
shear**

$$\gamma_t(R) = \Delta\Sigma(R)\,
\langle\Sigma_{\rm crit}^{-1}\rangle,$$

with the one-halo and projection contributions combined **additively**
(see the composition in {doc}`../science/index`).

The earlier **reduced-shear** formulation,

$$g_t(R) = \frac{\gamma_t(R)}{1 - \kappa(R)},$$

was **retired** (2026-05-11): its denominator makes the observable
nonlinear in $\Sigma$, which prevents the simple additive decomposition of
the one-halo and projection terms that the current pipeline relies on.
Restoring it would require composing $\Sigma(R)$ and $\Delta\Sigma(R)$
across all contributions *before* forming the ratio — i.e. a per-cluster
(or at least per-bin, population-averaged)
$\kappa(R) = \Sigma(R)\,\Sigma_{\rm crit}^{-1}$ assembled from the same
components. This is documented as a historical formulation, not a current
or planned output.

This section also fixes the distinction between a **per-cluster profile**
(e.g. $\langle\gamma_t^{1h}\rangle = N_{ij}[\gamma_t^{1h}]/N_{ij}[1]$) and
quantities already integrated over the selected population
(the raw $N_{ij}[\cdot]$ operators).

## Survey geometry

The effective survey area $\Omega(z)$ enters the number-count integrals
but is deliberately **excluded** from surface-density observables
(it cancels between numerator and normalisation for $\Sigma$-type
quantities — confirmed empirically: including it breaks the fiducial
likelihood closure).

```{todo}
Document: the $\Omega(z)$ model and its data files; mask treatment;
redshift dependence; any richness dependence; and the survey-area
convention to be adopted consistently on the Python side
(`richness_selection.survey_area` SurveyArea dataclass) when the mock
data vector is next regenerated.
```
