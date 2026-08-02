# Exact-gradient local-signature adversary

This experiment differentiates the complete local-triad signature objective
analytically. It replaces random screening with Riemannian gradient ascent on
the unit-energy, divergence-free, real Fourier-state sphere.

## Gradient certificate

For signature transfers `V_sigma`, define

```text
S = sum_sigma V_sigma,
R = sum_sigma |V_sigma|^2,
A_sig = |S| / sqrt(R).
```

Away from `S=0` and `R=0`, the exact signature weight in the derivative is

```text
dA_sig/dV_sigma = sign(S)/sqrt(R) - |S| V_sigma/R^(3/2).
```

`LocalSignatureObjective` propagates these weights through every complex
Fourier interaction, including the advecting, advected, and target slots.
The integrated self-test reports

```text
signature-ledger reconstruction error:       5.14e-19
A_sig gradient central-difference error:     6.75e-13
absolute-transfer gradient difference error: 7.61e-13
```

## Pointwise amplification obstruction

Twelve parallel flat-spectrum restarts, ten exact-gradient steps per restart,
give

```text
K       2        3         4         5         6
A_sig  4.5917  10.0918   19.6534   32.0093   49.3891
```

The fitted K2--K6 exponent is `2.16843`, far above the required `1/2`.
Therefore the pointwise LSF-4 amplification hypothesis is rejected. The
earlier random search missed these coherent states.

The machine artifact is
[`../local-signature-gradient-flat-K2-K6.json`](../local-signature-gradient-flat-K2-K6.json).

## Coupled magnitude test

Large amplification is coupled to a small signature-transfer square sum. On
the amplification-maximizing states, the normalized LSF-2 square-sum factor
falls from `0.004895` at K2 to `0.0004476` at K6. The product, however, must be
tested directly.

A separate exact-gradient run maximizes `|S|` itself. After forty accepted
steps per best restart it reports

```text
K             2        3        4        5         6
|S|         3.0785  11.5651  32.0015  77.7161  162.751
|S|/K^4     0.1924   0.1428   0.1250   0.1243    0.1256
```

The best runs still accept their fortieth step. The apparent K4 plateau is
critical, not closing, and no subcritical pointwise claim is made. Exact
fixed-energy concentration scaling already requires the unrestricted
asymptotic search to allow the classical K^(9/2) obstruction.

The transfer artifact is
[`../local-signature-gradient-transfer-flat-K2-K6.json`](../local-signature-gradient-transfer-flat-K2-K6.json).

## Dynamic restart point

The exact identity that remains useful is

```text
S^4/(Z P^3) = A_sig^4 R^2/(Z P^3).
```

On the replayable broad K3 trajectory, evaluated on the identical state
through cutoffs K3--K6, the critical integral is `4.39122e-6`, the bare
square-signature integral is `3.08593e-6`, and the last cutoff difference is
`1.51e-15`. The factorization residual is below `3.7e-19`. Repeating the same
time horizon with half the RK4 step changes the critical integral by only
`4.44e-7` relatively.

Thus neither factor admits the desired pointwise estimate in isolation. The
remaining candidate is a trajectory-integrated estimate for their product,
using viscosity and a fixed smooth initial datum. Its artifact is
[`../local-signature-trajectory-broad-K3-K6.json`](../local-signature-trajectory-broad-K3-K6.json).

Reproduce the two optimizations with

```bash
./build/navier_stokes_lab local-signature-gradient \
  --min-cutoff 2 --max-cutoff 6 --profile flat \
  --objective amplification --restarts 12 --workers 12 \
  --iterations 10 --line-search 16 --step 0.1 \
  --certificate proof/l4/analysis/local-signature-gradient/local-signature-gradient-flat-K2-K6.json

./build/navier_stokes_lab local-signature-gradient \
  --min-cutoff 2 --max-cutoff 6 --profile flat \
  --objective transfer --restarts 12 --workers 12 \
  --iterations 40 --line-search 18 --step 0.1 \
  --certificate proof/l4/analysis/local-signature-gradient/local-signature-gradient-transfer-flat-K2-K6.json
```
