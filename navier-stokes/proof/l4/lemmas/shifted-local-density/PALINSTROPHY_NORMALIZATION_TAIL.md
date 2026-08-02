# Palinstrophy-normalization tail tradeoff

This note isolates the smallest currently visible palinstrophy-normalization
subproblem in RQ-11.  The algebraic factorization and the zero-padding lemma
below are proved exactly.  The cutoff-uniform tail estimate is not proved, so
this note does not prove RQ-11, the local SLD lemma, L4 regularity, or the Clay
problem.

## Exact core--tail factorization

Fix a primitive-projective core of maximum height `H`.  Write the selected
bilinear field as

```text
b = b_c + b_t,
s_c = <Au,b_c>,       s_t = <Au,b_t>,
t_c = <Au,A b_c>,     t_t = <Au,A b_t>.
```

The part of the palinstrophy normalization outside the core is exactly

```text
(s_c+s_t)(t_c+t_t)-s_c t_c
    = s_c t_t+s_t t_c+s_t t_t.                    (PNT-1)
```

After multiplication by the power-one scale, its absolute objective is

```text
N_H(u)=3 |S_full| |s_c t_t+s_t t_c+s_t t_t|
       / (2 Z^2 P^3).                              (PNT-2)
```

The engine evaluates all three products independently and reconstructs
PNT-1.  The maximum relative reconstruction error in the current K8 and K12
certificates is below `3e-19`.

Let

```text
C1=||A^(1/2)b_c||_2^2,  C2=||A b_c||_2^2,
T1=||A^(1/2)b_t||_2^2,  T2=||A b_t||_2^2.
```

Cauchy--Schwarz gives the exact finite upper weights

```text
|s_c t_t| <= sqrt(Z C1 P T2),
|s_t t_c| <= sqrt(Z T1 P C2),
|s_t t_t| <= sqrt(Z T1 P T2).                      (PNT-3)
```

For every nonzero factor, actual/bound is exactly the product of the
corresponding H1 stretching alignment and H2 palinstrophy alignment.  This is
an identity, not a decay theorem.  Directly maximizing the tail stretching
alignment produced a highly aligned but very rough state with a much smaller
PNT-2 value, so an alignment-only lemma is false as a proof route.

## Exact zero-padding lemma

Let `u_K` be supported in a Galerkin cube of cutoff `K`, and embed it by zero
padding into any larger cube `K'`.  Then `Z`, `P`, `S_full`, every old-target
Fourier coefficient of `b_c` and `b_t`, and all four pairings `s_c,s_t,t_c,t_t`
are unchanged.  New bilinear output modes can occur, but `Au_K` is zero on
those modes, so they contribute zero to the pairings.  Consequently every
term in PNT-1 and PNT-2 is exactly invariant under zero padding.

The full norms `T1,T2` in PNT-3 need not be invariant: new output modes can
increase them without changing the pairings.  Therefore growth of the raw
Cauchy weight under zero padding is not evidence that PNT-2 grows.  The C++
self-test checks the open objective and all three factors under an
energy-preserving cutoff lift.

## Current stress result

Exact component gradients are alternated with the exact gradient of PNT-2.
At each cycle the `dominant` strategy reevaluates all three components and
selects the largest one before returning to the open-sum objective.

```text
state                                      PNT-2
original H8 K8 branch                      0.000345644731
K8 after ten alternating cycles            0.000573210985
same state zero-padded to complete K12      0.000573210985
K12 after one component/open gradient step  0.000575812137
```

At the tenth-cycle K8 state, the three PNT-1 contributions are

```text
core stretching x tail cross   0.0000941045   16.4%
tail stretching x core cross   0.0000457331    8.0%
tail stretching x tail cross   0.0004333734   75.6%
```

The joint Cauchy bound is `0.00844764` at K8.  Zero padding to K12 leaves the
actual value unchanged but raises that bound to `0.00849852`; actual/bound
falls from `0.06785` to `0.06745`.  One full K12 step improves the actual
objective by only `0.454%`.  These are finite stress tests, not upper bounds.

## Narrow candidate lemma

Every fixed projective ray is already controlled by the finite-family
plane--sphere argument.  For any fixed finite core, those estimates handle
the core--tail terms after summing finitely many core shapes.  The remaining
palinstrophy-normalization target is therefore the tail--tail estimate

```text
sup_(K,u) 3 |S_full| |<Au,b_t><Au,A b_t>|
           / (2 Z^2 P^3) < infinity.               (PNT-4)
```

PNT-4 must be proved jointly.  Bounding `T2` alone asks for an unavailable
higher derivative, while bounding the alignments alone was falsified by the
alignment adversary.  A successful proof must retain the tradeoff among the
H1/H2 tail size, the two pairings, and the common `Z,P,S_full` normalization.
Even PNT-4 would close only the open palinstrophy-normalization part; the
global response weight and enstrophy normalization in RQ-11 still require
their stated estimates.

## Reproduction

```bash
./build/navier_stokes_lab local-sld-projective-normalization-alternating \
  --state proof/l4/states/local-projective-normalization-alternating/H8-K8-tail-tail-eighth-cycle/K8.tsv \
  --output-state /tmp/K8-dominant.tsv \
  --certificate /tmp/K8-dominant.json \
  --selection double-triple-remainder-without-123 \
  --component dominant --projective-core-height 8 --cycles 2 \
  --component-iterations 12 --open-iterations 16 \
  --line-search 12 --lbfgs-history 8 --threads 12 --step 0.1

./build/navier_stokes_lab local-sld-projective-height-matrix \
  --state /tmp/K8-dominant.tsv --threads 12 \
  --exclude-triple-family --exclude-123 \
  --certificate /tmp/K8-dominant-matrix.json
```

Primary certificates:

- [`../../adversary/shifted-local-density/projective-normalization-alternating/H8-K8-dominant-cycles-nine-ten.json`](../../adversary/shifted-local-density/projective-normalization-alternating/H8-K8-dominant-cycles-nine-ten.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K8-open-palinstrophy-H8-dominant-tenth-cycle-matrix.json`](../../analysis/shifted-local-density/remainder-quartet/K8-open-palinstrophy-H8-dominant-tenth-cycle-matrix.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K12-zero-pad-open-palinstrophy-H8-dominant-tenth-cycle-matrix.json`](../../analysis/shifted-local-density/remainder-quartet/K12-zero-pad-open-palinstrophy-H8-dominant-tenth-cycle-matrix.json)
- [`../../adversary/shifted-local-density/projective-normalization-alternating/H8-K12-dominant-one-step.json`](../../adversary/shifted-local-density/projective-normalization-alternating/H8-K12-dominant-one-step.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K12-open-palinstrophy-H8-dominant-one-step-matrix.json`](../../analysis/shifted-local-density/remainder-quartet/K12-open-palinstrophy-H8-dominant-one-step-matrix.json)
