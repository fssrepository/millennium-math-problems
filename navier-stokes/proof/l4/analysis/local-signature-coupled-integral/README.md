# Coupled local-signature critical integral

This experiment attacks the remaining local L4 candidate directly. It does
not bound signature amplification and square-summed transfer separately. The
optimized quantity is

```text
J_local,N(u0) = integral_0^T |V_local,N(t)|^4 /
                              (Z_N(t) P_N(t)^3) dt.
```

Equivalently, the local squared-length signature ledger gives the exact
factorization

```text
|V_local|^4/(Z P^3)
    = A_sig^4 (sum_sigma |V_sigma|^2)^2/(Z P^3).
```

The pointwise factors are already known to fail separately. This run tests
their coupled trajectory product with the exact checkpointed RK4 adjoint.

## Independent dynamic multistart

`DynamicAdversary` owns one Galerkin, trajectory, adjoint, and L-BFGS context.
`DynamicAdversaryEnsemble` runs independent contexts concurrently, avoiding
shared mutable numerical state. Half of the additional starts perturb the
smooth continuation and half are independent random states. Every start is
retracted to fixed energy and the configured initial Sobolev cap.

At `E(0)=1`, homogeneous `H4^2<=100`, `nu=0.1`, `T=0.01`, and `dt=0.001`,
twelve dynamic starts and eight L-BFGS iterations per cutoff give

```text
K                    3             4             5             6
J_local refined  4.392286e-6  4.392821e-6  4.392876e-6  4.392883e-6
dt relative err  4.496e-7     4.562e-7     4.566e-7     4.568e-7
top-shell E      2.475e-7     2.130e-8     3.557e-11    1.133e-13
```

All four winners are restart zero, the continued smooth branch. At K6 the
best independent competitor reaches `4.392753e-6`, only `2.91e-5` below the
winner. Thus independent starts return to the same basin rather than exposing
a cutoff-concentrating branch. The K5-to-K6 state projection residual is
`7.12e-5`, and the initial homogeneous H4 squared norm decreases to `7.20317`.

The full run used an average `9.61` CPU cores, peaked at `345 MB` resident
memory, and took `7m04s`. It is finite adversarial evidence, not a uniform
analytical estimate and not a proof of global regularity.

Artifacts:

- [`../../adversary/local-signature-coupled-integral-h4-cap100-K3-K6-multistart12-lbfgs8.json`](../../adversary/local-signature-coupled-integral-h4-cap100-K3-K6-multistart12-lbfgs8.json)
- [`multistart-family.json`](multistart-family.json)

Reproduce with

```bash
./build/navier_stokes_lab adversary \
  --cutoffs 3,4,5,6 --restarts 12 --dynamic-restarts 12 \
  --generations 1 --dynamic-generations 8 --mutation 0.05 \
  --dynamic-objective critical-local-integral \
  --dynamic-optimizer gradient --gradient-method lbfgs \
  --sobolev-order 4 --sobolev-cap 100 \
  --nu 0.1 --evolve-time 0.01 --dt 0.001 \
  --threads 12 --backend fft \
  --dynamic-warm-state \
    proof/l4/states/helical-heterochiral-broad-spread/dynamic/K3.tsv
```

## Logical restart point

The observed branch is smooth and cutoff-stable, so it does not falsify the
trajectory-integrated coupled-product lemma. The next computational axes are
longer horizons and lower viscosity. The proof still requires a conventional
cutoff-uniform bound depending only on the fixed smooth initial datum,
viscosity, and finite time, without assuming a future high-Sobolev norm.

## Horizon continuation

Starting from the K6 multistart winner, a second twelve-start run at `T=0.02`
and six L-BFGS iterations gives

```text
J_local(0.02) / (2 J_local(0.01)) = 0.9974001829
log(Q(0.02)/Q(0)) /
    (2 log(Q(0.01)/Q(0)))       = 0.9655508952
dt relative error               = 4.5068e-7
```

The critical integral therefore accumulates slightly slower than linearly on
this finite branch. Restart zero wins again; the strongest independent start
is `8.09e-4` below it relatively. The artifact is
[`../../adversary/local-signature-coupled-integral-h4-cap100-K6-T002-multistart12-lbfgs6.json`](../../adversary/local-signature-coupled-integral-h4-cap100-K6-T002-multistart12-lbfgs6.json).

## Viscosity continuation

Reducing viscosity from `0.1` to `0.02` at K6 and `T=0.02`, with twelve new
dynamic starts and four L-BFGS iterations, gives

```text
J_local(nu=0.02) / J_local(nu=0.1) = 1.0038753823
log-Q gain ratio                       = 1.3093501057
dt relative error                     = 4.9732e-7
```

The critical integral rises by only `0.39%` on this branch when dissipation is
reduced fivefold, although Q grows more rapidly. The strongest independent
start is `2.43e-5` below the continued winner. This remains a short-time,
finite-cutoff observation. The artifact is
[`../../adversary/local-signature-coupled-integral-h4-cap100-K6-T002-nu002-multistart12-lbfgs4.json`](../../adversary/local-signature-coupled-integral-h4-cap100-K6-T002-nu002-multistart12-lbfgs4.json).

## Factor-correlation mechanism

On the four separately optimized K3--K6 trajectories,

```text
corr(log(A_sig^4), log(R^2/(Z P^3)))
  = -0.999980, -0.999978, -0.999978, -0.999978.
```

The amplification peaks at the final sample, while the square-signature and
coupled critical densities peak initially. At K6, `nu=0.02`, and `T=0.02`,
the correlation is `-0.999848`. The exact factorization residual remains
below `3.7e-19`. These branchwise observations explain why maximizing either
pointwise factor alone was misleading.

The corresponding universal derivative lemma is false. The
`local-signature-factor` command probed 240 independent decaying, flat, and
outer-half-flat states through K6. Forty-one samples have simultaneous
positive growth of `A_sig^4` and `R^2/(ZP^3)`. The largest common log-growth
rate is `19.8958`, and the largest critical-density log-growth rate is
`4391.14`. Halving the RK4 probe step preserves every simultaneous-growth
count and gives `19.8964` for the former rate. Therefore opposite factor
motion is not an algebraic Navier--Stokes identity; it is structure selected
by the current smooth extremizing branch.

Artifacts:

- [`nu002-K6-factor-correlation.json`](nu002-K6-factor-correlation.json)
- [`factor-adversary-K2-K6.json`](factor-adversary-K2-K6.json)
- [`factor-adversary-K2-K6-dt5e-5.json`](factor-adversary-K2-K6-dt5e-5.json)

Reproduce the falsification with

```bash
./build/navier_stokes_lab local-signature-factor \
  --min-cutoff 2 --max-cutoff 6 --samples 16 --workers 12 \
  --nu 0.1 --dt 0.0001 \
  --certificate proof/l4/analysis/local-signature-coupled-integral/factor-adversary-K2-K6.json
```
