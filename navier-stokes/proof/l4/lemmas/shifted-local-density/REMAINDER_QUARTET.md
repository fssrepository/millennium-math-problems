# Closed signature-remainder target

This note isolates the next open block after the complete `(m,m,2m)` family.
It gives the exact frequency loss of a dense energy-only estimate and the
precise collective gain required to close the remainder. It does not prove
the remainder bound or the full local SLD lemma.

## Fixed signatures already close

For a fixed squared-length signature `(r,s,t)`, fixing an input or target
leaves a plane--sphere intersection with degree `O(R)`. Therefore

```text
||B_sigma(v,w)||_2 <= C R^(3/2)||v||_2||w||_2.
```

The same quartet power count as in the doubling proof gives `R^5 E^2`, a
half derivative below the LQC-3 shell target `R^(11/2)E^2`. Consequently
every fixed remainder signature has a cutoff-independent estimate. The open
problem is summing the growing signature family, not any individual triple.

## Dense-count obstruction

For the complete thick local shell, a fixed target can interact with `O(R^3)`
inputs. Target-wise Cauchy--Schwarz then gives only

```text
||B_rem(u,u)||_2 <= C R^(5/2) E.                     (RQ-1)
```

The structural quartet entries have scale `R^7E^2`. The normalization entries
have the same power because

```text
S_rem = <A u,B_rem>             scales as R^(9/2)E^(3/2),
T_rem = <A B_rem,A u>           scales as R^(13/2)E^(3/2),
S_rem^2/Z and S_rem T_rem/P     scale as R^7E^2.     (RQ-2)
```

Thus a dense count loses exactly

```text
R^7/R^(11/2) = R^(3/2).                              (RQ-3)
```

This does not disprove the remainder lemma. It proves that energy, locality,
and an unsigned dense interaction count cannot establish it.

## Exact collective target

If an effective incidence estimate has degree `R^d`, then the bilinear
frequency power is `1+d/2`, and its structural quartet power is

```text
2 + 2(1+d/2) = 4+d.
```

Matching the LQC-3 exponent `11/2` requires

```text
d <= 3/2.                                             (RQ-4)
```

The raw local degree is `d=3`, so the remainder needs an `R^(3/2)` reduction
from signed signature coupling, an equivalent square-function inequality, or
a trajectory-propagated frequency envelope. Merely estimating every
signature and multiplying by their count reintroduces the full loss.

`RemainderQuartetClosure` evaluates all powers as exact rationals and the
self-test requires dense bilinear power `5/2`, dense quartet power `7`, target
power `11/2`, loss `3/2`, and required effective degree `3/2`.

## Current falsification screens

The two-scale axis/response family with `E_high=L^(-11/4)` does not create
growth in the remainder target ratio. In a 143-state full block scan the
largest absolute remainder ratio is `0.022443289746` at `L=3`; the per-scale
maximum decreases to `0.009228029403` at `L=12`. The mixed block vanishes on
this exactly self-similar family. These facts rule out simple dilation as the
missing obstruction, but they are finite evidence only.

An exact-gradient `lqc3-ratio` objective now optimizes

```text
|K_rem+G_rem|^2/[Z^(5/2)P^(3/2)]
```

directly. It is distinct from the older static `LQC-7` normalization and is
used to search for the collective dense-signature mechanism in RQ-4.

The first 12-start continuation gives

```text
K                         1          2          3          4
|K_rem+G_rem|/target   0.185413   0.194586   0.194707   0.194728
```

The squared-objective fitted slope over K2--K4 is `0.00217`. The K4 projected
gradient is still nonzero, so this is not a global optimum, but the warm branch
shows rapid cutoff flattening rather than the dense-count `R^(3/2)` growth.
The saved states make the result replayable.

## Reproduction

```bash
./build/navier_stokes_lab remainder-quartet-certificate \
  --certificate proof/l4/analysis/shifted-local-density/remainder-quartet/scaling-obstruction.json

./build/navier_stokes_lab local-sld-doubling-scale-scan \
  --min-scale 2 --max-scale 12 \
  --angle-min -1.2 --angle-max 1.2 --angle-count 13 \
  --energy-decay-power 2.75 --threads 12 \
  --certificate proof/l4/analysis/shifted-local-density/doubling-quartet/two-scale-full-block-L2-L12-angle-scan.json

./build/navier_stokes_lab local-closure-adversary \
  --objective lqc3-ratio --selection doubling-remainder \
  --min-cutoff 1 --max-cutoff 4 --restarts 12 --workers 12 \
  --iterations 12 --method lbfgs --backend direct \
  --certificate proof/l4/adversary/shifted-local-density/remainder-lqc3/K1-K4.json \
  --state-dir proof/l4/states/local-lqc3-ratio/doubling-remainder-K1-K4
```

Artifacts:

- [`../../analysis/shifted-local-density/remainder-quartet/scaling-obstruction.json`](../../analysis/shifted-local-density/remainder-quartet/scaling-obstruction.json)
- [`../../analysis/shifted-local-density/doubling-quartet/two-scale-full-block-L2-L12-angle-scan.json`](../../analysis/shifted-local-density/doubling-quartet/two-scale-full-block-L2-L12-angle-scan.json)
- [`../../adversary/shifted-local-density/remainder-lqc3/K1-K4.json`](../../adversary/shifted-local-density/remainder-lqc3/K1-K4.json)

The exact remaining statement is RQ-4 or an analytically equivalent signed
quartet estimate. No finite scan is treated as that proof.
