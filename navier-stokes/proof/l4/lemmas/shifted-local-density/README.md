# Shifted local critical-density candidate

This note records a scale-compatible replacement for the rejected pointwise
monotonicity and unshifted multiplicative-growth routes.

Define

```text
C_N(t) = |V_N,local(t)|^4/(Z_N(t) P_N(t)^3),
B_0    = E_N(0) P_N(0),
k_0    = sqrt(Z_N(0)/E_N(0)).
```

The exact candidate is

```text
L4.2-S:
d/dt log(C_N(t)+B_0) <= A(u_0,nu,T) k_0 Z_N(t)              (SLD-1)
```

with `A` independent of the Galerkin cutoff. SLD-1 is open.

## Exact compatibility checks

Under velocity-amplitude scaling, both `C_N` and `B_0` have degree four.
Under Navier--Stokes spatial scaling their exponent is two:

```text
C_N: 4*3 - 1 - 3*3 = 2,
B_0: -1 + 3 = 2.
```

The logarithmic time derivative has scaling exponent two. Since
`k_0` has exponent one and `Z_N` has exponent one, the right-hand side of
SLD-1 also has exponent two. `ShiftedCriticalDensityLemma` verifies these
relations with exact rational arithmetic in every self-test.

## Conditional closure

If SLD-1 is proved, the Galerkin energy identity gives

```text
integral_0^T Z_N(t) dt <= E_N(0)/(2 nu),

C_N(t)+B_0
  <= (C_N(0)+B_0) exp(A k_0 E_N(0)/(2 nu)).
```

Consequently `sup_N integral_0^T C_N(t) dt` is finite for each fixed smooth
initial datum, provided every constant in SLD-1 is cutoff independent. This
would close the remaining local L4.2 block when combined with the proved
moving far-tail estimate and a compatible transition-band estimate.

## Numerical screen

The exact discrete adjoint optimizes

```text
log((C_N(T)+B)/(C_N(0)+B))
```

for a configured finite `B`. At `E(0)=1`, `H4^2<=100`, `nu=0.1`, and
`T=0.001`, the K6 shift sweep gives

```text
B             1e-5       1e-4       1e-3
gain         0.011141    0.004508    0.000938
gain/(T k0 Z) 10.657       4.135       0.826
C_N(0)       1.87e-5     9.74e-5     1.89e-4
Delta C_N    3.21e-7     8.92e-7     1.12e-6
```

Time-step errors are below `2.7e-14`. The dependence on `B` is substantial,
so a free numerical shift cannot be promoted to a lemma. `B_0=E(0)P(0)` is
the current analytical choice because it has the exact required homogeneity
and depends only on the fixed smooth initial datum. The next proof task is to
differentiate `C_N` analytically and decide whether its local triad terms can
be bounded by the right-hand side of SLD-1 without a future high-Sobolev norm.

## Exact state-dependent shift adversary

The optimizer also differentiates `B_0=E(0)P(0)` as part of the initial state;
it is not frozen as a numerical parameter. The extra gradient term is

```text
dB_0 = P(0) dE(0) + E(0) dP(0).
```

The complete discrete gradient agrees with centered differences to
`2.06e-10`. Twelve starts and eight L-BFGS iterations give

```text
K                                  3          4          5          6
log((C(T)+B0)/(C(0)+B0))      7.918e-7   8.310e-7   8.558e-7   8.632e-7
gain/(T k0 Z(0))              7.240e-4   7.533e-4   7.751e-4   7.810e-4
```

Because different restarts win these four searches, the K6 winner was then
replayed independently at every lower cutoff. Finally, the identical K6 state
was zero-padded to K7 and K8 without optimization:

```text
K                                  6             7             8
shifted log gain              8.632216851e-7 8.632152689e-7 8.632152671e-7
gain/(T k0 Z(0))              7.810352793e-4 7.810294739e-4 7.810294723e-4
relative change                    -        -7.43e-6      -2.14e-9
```

This finite branch is time-step stable and cutoff converged. It does not prove
SLD-1 or bound its global constant, but it does not falsify the candidate.

## Instantaneous derivative oracle

`shifted-density` evaluates the left side of SLD-1 directly, without an RK4
horizon:

```text
dC_local/dt = <gradient C_local, Navier--Stokes RHS>,
R_SLD       = (d/dt log(C_local+E0P0))/(k0 Z(0)).
```

It evaluates the K6 winner in `0.22 s`. Direct-triad and FFT RHS backends give
identical serialized values. The separately optimized winners give

```text
K                  3             4             5             6
R_SLD       7.239567e-4  7.534130e-4  7.754572e-4  7.815121e-4
```

At K6 the instantaneous rate differs from the `T=0.001` average by only
`6.10e-4` relatively. This command is now the fast falsification path for
SLD-1 states; trajectory integration is reserved for candidates that survive
the instantaneous screen.

Reproduce with

```bash
./build/navier_stokes_lab shifted-density \
  --state proof/l4/states/local-critical-ep-log-gain-multistart/dynamic/K6.tsv \
  --nu 0.1 --threads 12 --backend fft \
  --certificate /tmp/K6-shifted-density.json
```

Artifacts:

- [`../../adversary/local-critical-shifted-log-gain1e-5-h4-cap100-K6-multistart6-lbfgs8.json`](../../adversary/local-critical-shifted-log-gain1e-5-h4-cap100-K6-multistart6-lbfgs8.json)
- [`../../adversary/local-critical-shifted-log-gain1e-4-h4-cap100-K3-K6-multistart12-lbfgs8.json`](../../adversary/local-critical-shifted-log-gain1e-4-h4-cap100-K3-K6-multistart12-lbfgs8.json)
- [`../../adversary/local-critical-shifted-log-gain1e-3-h4-cap100-K6-multistart6-lbfgs8.json`](../../adversary/local-critical-shifted-log-gain1e-3-h4-cap100-K6-multistart6-lbfgs8.json)
- [`../../adversary/local-critical-ep-log-gain-h4-cap100-K3-K6-multistart12-lbfgs8.json`](../../adversary/local-critical-ep-log-gain-h4-cap100-K3-K6-multistart12-lbfgs8.json)
- [`../../adversary/local-critical-ep-log-gain-K6-projections-K3-K6.json`](../../adversary/local-critical-ep-log-gain-K6-projections-K3-K6.json)
- [`../../adversary/local-critical-ep-log-gain-K6-lifts-K6-K8.json`](../../adversary/local-critical-ep-log-gain-K6-lifts-K6-K8.json)
- [`../../analysis/local-critical-increase/K6-instantaneous-E0P0-shifted-density.json`](../../analysis/local-critical-increase/K6-instantaneous-E0P0-shifted-density.json)
