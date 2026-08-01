# Local critical-density increase

This experiment tests whether the remaining coupled local density could be
controlled by a pointwise monotonicity lemma. Define

```text
C_local(u) = |V_local(u)|^4/(Z(u) P(u)^3).
```

The optimized objective is the exact endpoint difference

```text
Delta C_local = C_local(u(T)) - C_local(u(0)).
```

`SpectralAdjoint::critical_increase_gradient` propagates the terminal gradient
through the discrete RK4 map and subtracts the initial gradient. The
deterministic directional test agrees with centered finite differences to
`5.23e-12` relatively.

## K3--K6 result

Parameters: `E(0)=1`, homogeneous `H4^2<=100`, `nu=0.1`, `T=0.001`,
`dt=0.0005`, twelve independent dynamic starts, eight projected L-BFGS
iterations, and twelve CPU workers.

```text
K                     3             4             5             6
C_local(0)       2.228636e-4  2.261852e-4  2.213693e-4  2.234314e-4
C_local(T)       2.238284e-4  2.272846e-4  2.224945e-4  2.245664e-4
Delta C_local    9.648047e-7  1.099410e-6  1.125187e-6  1.134951e-6
top-shell E      5.013e-4     1.253e-4     2.259e-5     1.408e-6
```

Every increase is positive. Repeating each winning trajectory with half the
time step changes the objective by less than `2.5e-14` relatively. The fitted
top-shell energy exponent is `-8.18`; the K5-to-K6 projective residual is
`1.74e-3`.

Therefore the universal statement "`C_local` is nonincreasing" is false
(F010). The branch is smooth and its increase flattens with cutoff, so the
calculation does not show a singularity and does not falsify the desired
time-integrated bound.

## K6 optimizer refinement

Starting from the eight-step K6 winner, a fresh twelve-start run with 24
L-BFGS iterations raises `Delta C_local` from `1.13495104199e-6` to
`1.13820060675e-6`, a relative gain of only `0.00286318`. The strongest
alternate restart is `7.91808e-4` lower than the continued winner. The
refined trajectory begins at `C_local=2.28456361874e-4`, remains on the
`H4^2=100` boundary, and passes time-step refinement at `2.68e-14` relative
error. The run averaged `9.08` CPU cores, used `281 MB` peak resident memory,
and completed in `3m30s`.

This separates most remaining optimizer progress from cutoff growth: much
longer optimization changes the objective by less than three tenths of one
percent.

## Horizon continuation and factor motion

At K6, doubling the horizon to `T=0.002` and running twelve starts with eight
new L-BFGS iterations gives

```text
Delta C_local(0.002)/(2 Delta C_local(0.001)) = 0.999368654
Delta C_local(0.002)/0.002                    = 0.00113748201
time-step relative error                      = 2.75e-14
```

The short-time increase is therefore almost exactly linear rather than
accelerating. A full signature replay of the winning trajectory gives

```text
corr(log(A_sig^4), log(R^2/(ZP^3))) = +0.999994857.
```

Both factors and the coupled density peak at the final RK4 sample, with exact
factorization residual `2.32e-19`. This is qualitatively different from the
integral-maximizing branch, whose two factors are almost perfectly
anticorrelated. The endpoint adjoint has therefore found a relevant
high-density simultaneous-growth branch, not merely the low-density random
counterexamples used for F009.

## Viscosity continuation

Repeating the K6, `T=0.002` search at `nu=0.02` gives

```text
Delta C_local(nu=0.02) / Delta C_local(nu=0.1) = 0.909393069
log-gain ratio                                 = 0.823480963
time-step relative error                       = 2.00e-15
```

The optimizer improves its lower-viscosity warm start by `2.48%`, but the
final absolute and relative growth remain smaller than at `nu=0.1`. On this
branch the transient is therefore not amplified by simply reducing
dissipation fivefold.

## Relative-growth objective and its zero-set obstruction

The exact adjoint also differentiates

```text
log(C_local(u(T))/C_local(u(0))).
```

Its centered-difference error is `2.08e-11`. An unconstrained K3--K6 search
does not return to the high-density branch. Instead it drives
`C_local(0)` to the numerical positivity floor near `1e-30`, while
`C_local(T)` is about `1e-17`, producing log gains from `28.0` to `31.81` and
projected gradients as large as `3.8e7`. K5 and K6 are already projectively
close, but their objective is limited by the denominator guard rather than a
physical cutoff mechanism.

This is a useful obstruction: a purely multiplicative Gronwall ansatz is
singular near the zero-transfer set and the raw log objective is not a valid
high-density proxy. The engine therefore supports
`--critical-density-shift B`, which instead optimizes

```text
log((C_local(T)+B)/(C_local(0)+B)).
```

The shifted discrete gradient passes centered differences at `2.79e-10` for
`B=1e-6`. A shift sweep can now test additive source terms without rewarding
numerically vanishing initial transfer.

For `B=1e-4`, the twelve-start K3--K6 continuation gives

```text
K                              3         4         5         6
shifted log gain          0.004269  0.004445  0.004495  0.004508
gain/(T k0 Z(0))              3.889     4.068     4.119     4.135
C_local(0)                 9.06e-5   9.73e-5   9.80e-5   9.74e-5
Delta C_local              8.15e-7   8.79e-7   8.92e-7   8.92e-7
```

Restart zero wins at every cutoff. Top-shell energy decays with fitted
exponent `-7.48`, and the K5-to-K6 projection residual is `4.14e-3`. This
shifted branch is again smooth and cutoff-stabilizing; it supplies a finite
candidate constant near `4.14`, not a proof that such a constant is universal.

The scale-compatible state-dependent choice `B_0=E(0)P(0)` is implemented by
a separate exact adjoint. Its normalized K3--K6 values are
`7.240e-4, 7.533e-4, 7.751e-4, 7.810e-4`. Replaying the K6 winner at K7 and K8
changes its objective by `-7.43e-6` and then `-2.14e-9` relatively. This branch
therefore survives the finite screen without cutoff growth. The mathematical
candidate and conditional closure are stated in
[`../../lemmas/shifted-local-density/README.md`](../../lemmas/shifted-local-density/README.md).

Artifacts:

- [`../../adversary/local-critical-increase-h4-cap100-K3-K6-multistart12-lbfgs8.json`](../../adversary/local-critical-increase-h4-cap100-K3-K6-multistart12-lbfgs8.json)
- [`../../adversary/local-critical-increase-h4-cap100-K6-multistart12-lbfgs24.json`](../../adversary/local-critical-increase-h4-cap100-K6-multistart12-lbfgs24.json)
- [`../../adversary/local-critical-increase-h4-cap100-K6-T002-multistart12-lbfgs8.json`](../../adversary/local-critical-increase-h4-cap100-K6-T002-multistart12-lbfgs8.json)
- [`../../adversary/local-critical-increase-h4-cap100-K6-T002-nu002-multistart12-lbfgs8.json`](../../adversary/local-critical-increase-h4-cap100-K6-T002-nu002-multistart12-lbfgs8.json)
- [`../../adversary/local-critical-log-gain-h4-cap100-K3-K6-multistart12-lbfgs8.json`](../../adversary/local-critical-log-gain-h4-cap100-K3-K6-multistart12-lbfgs8.json)
- [`../../adversary/local-critical-shifted-log-gain1e-4-h4-cap100-K3-K6-multistart12-lbfgs8.json`](../../adversary/local-critical-shifted-log-gain1e-4-h4-cap100-K3-K6-multistart12-lbfgs8.json)
- [`../local-signature-coupled-integral/critical-increase-family.json`](../local-signature-coupled-integral/critical-increase-family.json)
- [`log-gain-family.json`](log-gain-family.json)
- [`shifted-log-gain1e-4-family.json`](shifted-log-gain1e-4-family.json)
- [`K6-T002-factor-trajectory.json`](K6-T002-factor-trajectory.json)

Reproduce with

```bash
./build/navier_stokes_lab adversary \
  --cutoffs 3,4,5,6 --restarts 12 --dynamic-restarts 12 \
  --generations 1 --dynamic-generations 8 --mutation 0.05 \
  --dynamic-objective critical-local-increase \
  --dynamic-optimizer gradient --gradient-method lbfgs \
  --sobolev-order 4 --sobolev-cap 100 \
  --nu 0.1 --evolve-time 0.001 --dt 0.0005 \
  --threads 12 --backend fft \
  --dynamic-warm-state \
    proof/l4/states/helical-heterochiral-broad-spread/dynamic/K3.tsv
```
