# Local helical decomposition

`HelicalTriadLedger` decomposes every divergence-free Fourier coefficient into
curl eigenvectors

```text
i k cross h_s(k) = s |k| h_s(k),   s in {-1,+1},
```

and reconstructs vortex stretching from all eight sign sectors
`(s_p,s_q,s_k)`. The decomposition is exact: the current certificate has
maximum velocity, total-stretching, and local-stretching reconstruction
residuals below `2.5e-20`.

## Adversarial-state observation

For the K3 local-objective endpoint,

```text
V_local              = -2.2778266270e-1
V_local,homochiral   = -3.7146823897e-6
V_local,heterochiral = -2.2777894802e-1
E_plus / E_minus     = 0.4999997884 / 0.5000002116
```

Thus the signed local extremizer is almost entirely heterochiral. Its
homochiral absolute interaction sum is `6.46194e-1`, so the small signed value
comes from cancellation rather than absent homochiral interactions. The
separately optimized nonlocal endpoint has the same qualitative local split:
`V_local,homochiral=-5.18518e-4` and
`V_local,heterochiral=-3.41206e-2`.

## Rejected identity

The observed cancellation is not universal. Projecting deterministic random
divergence-free states onto pure positive or pure negative helicity leaves no
heterochiral sectors but produces nonzero local enstrophy stretching. The
certificate reaches

```text
max |V_local,homochiral| = 6.35947e-3.
```

Therefore candidate F006, exact homochiral cancellation, is rejected. The
surviving analytical direction is narrower: identify a quantitative dynamic
depletion of the heterochiral local sectors.

## Sector-targeted adversary

`HelicalSectorObjective` has an analytic gradient for the homochiral and
heterochiral local critical densities. Its central-difference errors are
`5.77e-12` for signed stretching and `6.71e-12` for the critical density.
Twelve parallel fixed-energy heterochiral restarts with eight accepted steps
raise the best static critical density from `6.04266e-9` to `1.33682e-4`.

`HelicalSectorAdjoint` differentiates the trapezoidal time integral through
every full Navier-Stokes RK4 step. Its trajectory-gradient error is
`3.56e-13`. The first 12-restart, four-step trajectory search raises the
heterochiral local critical integral from `1.91932e-11` to `5.07257e-8`.

These searches show that the heterochiral target is numerically active rather
than automatically depleted. The next cutoff continuation must determine
whether it approaches a smooth projective branch or grows with resolution.
The finite endpoint pattern is routing evidence, not a regularity proof.

## Cutoff continuation result

`HelicalAdversaryCli` now supports parallel multistarts and evaluation-only
runs. Starting from the converged K3 local state, eight trajectory-gradient
steps on each of twelve starts raised the original branch from
`4.38986521607e-6` to `4.39045255695e-6`; all perturbed starts were lower.
Lifting that exact state without further optimization and evaluating the same
`nu=0.1`, `dt=0.001`, `T=0.01` trajectory gives

```text
K:                       3                  4                  5                  6
J_heterochiral: 4.39045255695e-6   4.39045682754e-6   4.39045682741e-6   4.39045682741e-6
adjacent relative:       -          9.72697956e-7      2.98371299e-11     2.66978206e-15
```

The aggregate machine certificate is
`../../adversary/helical-heterochiral-cutoff-scan-K3-K6.json`. It shows that
this finite branch is a smooth low-mode branch, not a cutoff-concentration
branch.

A separate static search directly targets the fourth power of the required
depletion coefficient,

```text
|V_local,heterochiral|^4 / (Z P^3).
```

Sixteen accepted steps with twelve large (`0.5`) spectral perturbations at
each cutoff give `4.4026806e-4`, `4.4031124e-4`, and `4.4031558e-4` at K3, K4,
and K5. The unperturbed continuation won every multistart. This is useful
evidence that the current extremizer is numerically stable, but it is not a
global optimization certificate and cannot establish the required dynamic
bound.

The remaining analytical target is unchanged: prove a trajectory-restricted,
cutoff-independent bound for the heterochiral coefficient using only the
initial smooth datum, viscosity, and finite time. A universal static bound is
not a substitute for that dynamic statement.
