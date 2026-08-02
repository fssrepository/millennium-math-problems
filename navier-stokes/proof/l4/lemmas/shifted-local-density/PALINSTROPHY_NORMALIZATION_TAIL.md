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

With `s=s_c+s_t`, the same identity has the canonical two-term form

```text
s t-s_c t_c = s t_t+s_t t_c.                      (PNT-1b)
```

After multiplication by the power-one scale, its absolute objective is

```text
N_H(u)=3 |S_full| |s_c t_t+s_t t_c+s_t t_t|
       / (2 Z^2 P^3).                              (PNT-2)
```

The engine evaluates all three PNT-1 products and both PNT-1b products
independently.  The maximum relative reconstruction error in the current K8
and K12 certificates is below `4e-19`.

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

Writing `B1=||A^(1/2)b||_2^2`, the canonical high-derivative channel obeys

```text
|s t_t| <= sqrt(Z B1 P T2).                        (PNT-3b)
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
self-test checks the open objective and all four optimized factors under an
energy-preserving cutoff lift.

## Current stress result

Exact component gradients are alternated with the exact gradient of PNT-2.
At each cycle the `dominant` strategy reevaluates the two disjoint PNT-1b
components and selects the larger one before returning to the open-sum
objective.  The older three-way components remain available as diagnostics.

```text
state                                      PNT-2
original H8 K8 branch                      0.000345644731
K8 after ten alternating cycles            0.000573210985
same state zero-padded to complete K12      0.000573210985
K12 after one component/open gradient step  0.000575812137
K8 after one canonical two-term cycle        0.000590223320
K12 after one canonical component/open step  0.000592976913
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

At the eleventh-cycle K8 state, the disjoint PNT-1b values are

```text
selected stretching x tail cross  0.0005389546   91.3%
tail stretching x core cross      0.0000512687    8.7%
```

The two-term and three-term reconstruction errors are both `8.97e-20`.
Across the independently optimized H=8,16,32,64 scan, the selected/tail
channel is dominant in all four rows.  Its fitted finite height slope is
`-0.1498`, while its raw Cauchy-bound slope is `+0.0378`; neither slope is a
uniform theorem.
On the one-step K12 state the selected/tail channel is `0.000542046495`, or
`91.4%` of the open value; its two-term reconstruction error is `8.93e-20`.

## Exact alignment--majorant factorization

For the canonical selected/tail channel define

```text
A_H = s^2 t_t^2/(Z P B1 T2),
M_H = 3 |S_full| sqrt(Z B1 P T2)/(2 Z^2 P^3).
```

Then the selected-channel PNT value factors exactly as

```text
3 |S_full s t_t|/(2 Z^2 P^3) = M_H sqrt(A_H).       (PNT-6)
```

The engine now has independent exact-gradient objectives for both factors.
`projective-normalization-alignment-ratio` maximizes `A_H`, while
`projective-normalization-cauchy-ratio` maximizes

```text
M_H^2 = 9 S_full^2 B1 T2/(4 Z^3 P^5).              (PNT-7)
```

The PNT-2 K8 winner starts with `sqrt(A_H)=0.0832197` and
`M_H=0.00647629`.  Successive alignment optimization raises
`sqrt(A_H)` to `0.520422`, `0.689024`, and `0.771198`, but the selected PNT
value falls from `0.000538955` to `1.777e-7`, `1.084e-8`, and `1.205e-9`.
The last value is about 447,000 times below the selected-channel record.
A 65-sample normalized affine scan also decreases monotonically in PNT value
from the PNT endpoint toward the alignment endpoint.  These finite results
rule out treating near-Cauchy equality as a useful standalone proof target;
they do not prove a uniform tradeoff.

Direct majorant optimization gives the complementary stress test.  After 24
K8 L-BFGS steps, `M_H` rises from `0.00647629` to `0.0104101`, while
`sqrt(A_H)` falls to `0.0258487` and the selected PNT value becomes
`0.000269086`.  Exact zero padding of this state to complete K12 raises
`M_H` slightly further to `0.0104580`; one K12 gradient step reaches
`0.0104740`, with selected PNT value `0.000269574` and actual/bound ratio
`0.0257375`.  The K8 objective agrees bit-for-bit with the square of the
independent height-matrix Cauchy bound, and its analytic-gradient
central-difference error is `7.14e-12` in the self-test.  The cutoff behavior
of `M_H`, not alignment alone, is therefore the sharper computational target
for a PNT-5 lemma.

## Narrow candidate lemma

Every fixed projective ray is already controlled by the finite-family
plane--sphere argument, but that fact alone does not bound a fixed-core
factor paired with the complete tail.  The correct open target retains both
disjoint PNT-1b channels:

```text
sup_(K,u) 3 |S_full| |s t_t+s_t t_c|
           / (2 Z^2 P^3) < infinity.               (PNT-4)
```

The numerically dominant high-derivative subtarget is

```text
sup_(K,u) 3 |S_full s t_t|/(2 Z^2 P^3) < infinity. (PNT-5)
```

The complementary `S_full s_t t_c` estimate is still required and is not
marked proved.  Bounding `T2` alone asks for an unavailable higher derivative,
while bounding the alignments alone was falsified by the alignment adversary.
A successful proof must retain the tradeoff among the H1/H2 tail size, the
two pairings, and the common `Z,P,S_full` normalization.  Even PNT-4 would
close only the open palinstrophy-normalization part; the global response
weight and enstrophy normalization in RQ-11 still require their stated
estimates.

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

./build/navier_stokes_lab local-closure-adversary \
  --objective projective-normalization-cauchy-ratio \
  --selection double-triple-remainder-without-123 \
  --min-cutoff 8 --max-cutoff 8 --restarts 1 --workers 12 \
  --iterations 24 --method lbfgs --backend direct --lean \
  --preserve-warm-layout --projective-core-height 8 \
  --warm-state proof/l4/states/local-projective-normalization-alternating/H8-K8-two-term-eleventh-cycle/K8.tsv \
  --certificate /tmp/K8-normalization-cauchy.json \
  --state-dir /tmp/K8-normalization-cauchy
```

Primary certificates:

- [`../../adversary/shifted-local-density/projective-normalization-alternating/H8-K8-dominant-cycles-nine-ten.json`](../../adversary/shifted-local-density/projective-normalization-alternating/H8-K8-dominant-cycles-nine-ten.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K8-open-palinstrophy-H8-dominant-tenth-cycle-matrix.json`](../../analysis/shifted-local-density/remainder-quartet/K8-open-palinstrophy-H8-dominant-tenth-cycle-matrix.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K12-zero-pad-open-palinstrophy-H8-dominant-tenth-cycle-matrix.json`](../../analysis/shifted-local-density/remainder-quartet/K12-zero-pad-open-palinstrophy-H8-dominant-tenth-cycle-matrix.json)
- [`../../adversary/shifted-local-density/projective-normalization-alternating/H8-K12-dominant-one-step.json`](../../adversary/shifted-local-density/projective-normalization-alternating/H8-K12-dominant-one-step.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K12-open-palinstrophy-H8-dominant-one-step-matrix.json`](../../analysis/shifted-local-density/remainder-quartet/K12-open-palinstrophy-H8-dominant-one-step-matrix.json)
- [`../../adversary/shifted-local-density/projective-normalization-alternating/H8-K8-two-term-eleventh-cycle.json`](../../adversary/shifted-local-density/projective-normalization-alternating/H8-K8-two-term-eleventh-cycle.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K8-open-palinstrophy-H8-two-term-eleventh-cycle-matrix.json`](../../analysis/shifted-local-density/remainder-quartet/K8-open-palinstrophy-H8-two-term-eleventh-cycle-matrix.json)
- [`../../adversary/shifted-local-density/projective-normalization-alternating/H8-K12-two-term-one-step.json`](../../adversary/shifted-local-density/projective-normalization-alternating/H8-K12-two-term-one-step.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K12-open-palinstrophy-H8-two-term-one-step-matrix.json`](../../analysis/shifted-local-density/remainder-quartet/K12-open-palinstrophy-H8-two-term-one-step-matrix.json)
- [`../../adversary/shifted-local-density/projective-normalization-alignment/H8-K8-ninety-six-step-continuation.json`](../../adversary/shifted-local-density/projective-normalization-alignment/H8-K8-ninety-six-step-continuation.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K8-normalization-alignment-H8-ninety-six-step-continuation-matrix.json`](../../analysis/shifted-local-density/remainder-quartet/K8-normalization-alignment-H8-ninety-six-step-continuation-matrix.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K8-H8-normalization-joint-alignment-tradeoff-line.json`](../../analysis/shifted-local-density/remainder-quartet/K8-H8-normalization-joint-alignment-tradeoff-line.json)
- [`../../adversary/shifted-local-density/projective-normalization-cauchy/H8-K8-twenty-four-step.json`](../../adversary/shifted-local-density/projective-normalization-cauchy/H8-K8-twenty-four-step.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K8-normalization-cauchy-H8-twenty-four-step-matrix.json`](../../analysis/shifted-local-density/remainder-quartet/K8-normalization-cauchy-H8-twenty-four-step-matrix.json)
- [`../../adversary/shifted-local-density/projective-normalization-cauchy/H8-K12-zero-pad-evaluate.json`](../../adversary/shifted-local-density/projective-normalization-cauchy/H8-K12-zero-pad-evaluate.json)
- [`../../adversary/shifted-local-density/projective-normalization-cauchy/H8-K12-one-step.json`](../../adversary/shifted-local-density/projective-normalization-cauchy/H8-K12-one-step.json)
- [`../../adversary/shifted-local-density/projective-normalization-cauchy/H8-K12-one-step-selected-channel-evaluate.json`](../../adversary/shifted-local-density/projective-normalization-cauchy/H8-K12-one-step-selected-channel-evaluate.json)
