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

## Height-decay candidate and locality stress

The current explicit sufficient target for the selected/tail channel is

```text
sup_(K,u,H>=1) H^(1/2) M_H^2
 = sup_(K,u,H>=1) 9 H^(1/2) S_full^2 B1 T2/(4 Z^3 P^5)
 < infinity.                                               (PNT-8)
```

Equivalently, `M_H <= C H^(-1/4)`.  Since `0 <= A_H <= 1`, PNT-8 would
control the PNT-5 tail with explicit height decay.  It would not control the
second PNT-4 channel, and it is not proved.

Independent 24-step K8 optimizations at `H=8,16,32,64` give

```text
H                     8          16          32          64
M_H             0.0104101   0.00931419  0.00716035  0.00432474
H^(1/4) M_H     0.0175076   0.0186284   0.0170303   0.0122322
```

The fitted K8 height slope is `-0.41813`.  Exact zero padding to K12 changes
the fitted slope to `-0.33196`; one K12 gradient step at each height changes
it to `-0.27537`.  On those one-step K12 states, `H^(1/4)M_H` stays between
`0.01672` and `0.01899`, with compensated fitted slope `-0.02537`.  A separate
four-step H128 K12 stress run reaches `M_H=0.00396002` and
`H^(1/4)M_H=0.0133199`.  These values motivate PNT-8 but do not establish its
exponent or a uniform constant.

Zero padding preserves `Z`, `P`, and `S_full` exactly but increases the
aggregate norm factors.  At H=8,16,32,64, respectively, K8-to-K12 padding
raises `B1` by `0.16%,0.61%,2.84%,6.32%`, raises `T2` by
`0.76%,2.66%,12.36%,36.15%`, and raises `M_H` by
`0.46%,1.63%,7.49%,20.32%`.  Therefore a proof of PNT-8 cannot infer tail
decay from fixed-cutoff output support.

The most obvious low/high counterexample is blocked by the exact local-triad
condition

```text
max(|p|^2,|q|^2,|k|^2) <= 4 min(|p|^2,|q|^2,|k|^2).        (PNT-9)
```

A K2 low core plus one structured outer-shell satellite has `T2=M_H=0` for
H64 through K12: strongly separated low/high interactions never enter this
local remainder.  A dense comparable-frequency satellite does activate
high/high interactions.  Its H64 majorant is `1.32e-7,7.27e-7,2.13e-6,
2.76e-6,2.08e-6` at K=5,6,8,10,12.  The increase is a tail-entry threshold;
the K8--K12 rows plateau and then decrease, while actual/Cauchy is at most
`1.25e-5`.  Thus this bounded scan neither falsifies PNT-8 nor supplies the
missing proof.  It identifies the analytic bottleneck as a local high/high
band estimate for `T2`, coupled to the common `Z,P,S_full,B1` normalization.

## Exact tail Gram--Schur reduction

Let `b_j` be the aggregate bilinear output in the dyadic primitive-height
shell `j`, restricted to heights above H, and set

```text
d_j = ||A b_j||_2^2,
D_H = sum_(j>H) d_j,
R_H = max_i sum_j |<A b_i,A b_j>|/sqrt(d_i d_j).             (PNT-10)
```

Zero rows are omitted from the normalized sum.  The finite symmetric Schur
test gives the exact algebraic reduction

```text
T2 = ||A sum_(j>H)b_j||_2^2 <= R_H D_H.                    (PNT-11)
```

Consequently the joint estimate

```text
sup_(K,u,H) H^(1/2)
  9 S_full^2 B1 R_H D_H/(4 Z^3 P^5) < infinity             (PNT-12)
```

is sufficient for PNT-8.  A uniform bound on `R_H` together with the
corresponding diagonal estimate is a stronger sufficient split, not a proved
claim.

The new C++ ledger reconstructs `T2` both from the signed Gram sum and from
the independent Cauchy objective.  On all four independently optimized K8
states, every inter-shell H2 Gram contribution has the same sign, so
`T2/absolute-Gram=1`: cancellation is not responsible for the observed
height decay.  The maximum `R_H` is `1.65993`; the K12 one-step scan gives
`1.66034`.  The fitted squared-height slopes are

```text
scan       M_H^2       diagonal       Schur bound
K8        -0.83626      -0.65322       -0.88137
K12       -0.55074      -0.45692       -0.57874
```

At K12 the maximum compensated diagonal and Schur quantities are
`2.66683e-4` and `3.89149e-4`.  The finite Schur bound is 89--98% sharp on
these states.  This removes a possible cancellation explanation and narrows
the analytical problem to a number-theoretic/projective diagonal tail bound
plus control of `R_H`.  Neither finite maximum is a uniform theorem.

The pair ledger also exposed a more explicit but substantially stronger
sufficient route.  If the height-shell aggregates were to satisfy

```text
|<A b_i,A b_j>| <= C 2^(-|i-j|)
                    ||A b_i||_2 ||A b_j||_2,               (PNT-13)
```

then each normalized absolute Gram row would be bounded by
`1+2C sum_(g>=1)2^(-g)=1+2C`.  The original PNT-majorant extremizers looked
misleadingly favorable: their sharp finite
`C=max_(i!=j) 2^|i-j| |corr(i,j)|` values were `0.74069` at K8 and `0.74075`
at K12.

An exact-gradient adversary directly targeting the squared PNT-13 ratio
rejects this as a standalone proof mechanism.  It produces

```text
cutoff  shells  gap  raw correlation  2^gap |correlation|
K8       4,8     4      0.985710             15.7714
K8       2,8     6      0.998998             63.9359
K12      2,9     7      0.994971            127.3563
```

The weighted values are close to their algebraic maxima `2^gap`, not to the
earlier `0.74` plateau.  These finite adversaries invalidate PNT-13 as the
next proof route; a rigorous theorem that no uniform constant exists would
still require an explicit scalable family.

Crucially, the same states do not invalidate the necessary joint PNT-12
tradeoff.  Their maximum `H^(1/2)`-compensated Schur majorant is only
`3.56102e-13`, compared with `3.89149e-4` on the direct PNT-majorant stress
states.  The three normalization common factors
`9 S_full^2 B1/(4 Z^3 P^5)` are `1.03104e-21`, `1.59369e-17`, and
`2.19156e-24`; near-perfect shell alignment is bought by collapsing the
common normalized factor.  Thus PNT-12, not PNT-13, is the active lemma.

A mode-resolved ledger shows why the false route is easy to realize at finite
cutoff.  On the K8 gap-six state, output pairs `(8,-1,8)` and `(8,1,8)` carry
`73.798%` and `26.005%` of the signed Gram pairing.  On the K12 gap-seven
state, `(-12,12,0)` and `(12,0,12)` carry `85.731%` and `13.200%`.
The effective shared-output counts are only `1.63` and `1.33`.  The
correlation adversary therefore concentrates both shell outputs on one or two
Galerkin-wall modes; this is a structural clue, not an infinite-dimensional
bound.

Exact ordered-triad attribution shows that the low-dimensional output does
not come from an equally sparse high-shell input at increasing cutoff.  At
K8, the dominant `(8,-1,8)` output has 2 shell-2 and 714 shell-8 ordered
interactions, with effective counts `1.09` and `6.75`.  The shell-2 pair has
primitive shape `(1,2,3)`.  At K12, the dominant `(-12,12,0)` output has 19
shell-2 and 968 shell-9 interactions, with effective counts `4.11` and
`43.43`.  Thus the common output remains concentrated while the high-shell
generator becomes denser.  This rules out declaring a sparse scalable
counterexample from the present data and reinforces the need for the coupled
normalization in PNT-12.

## Direct coupled PNT-12 adversary

The engine now differentiates the actual fixed-row Schur quantity

```text
J_(H,i)(u) = H^(1/2)
  9 S_full^2 B1 R_(H,i) D_H/(4 Z^3 P^5),
R_(H,i) = sum_j |<A b_i,A b_j>|/sqrt(d_i d_j).       (PNT-14)
```

On a stratum where the finitely many Gram signs are fixed, this objective is
smooth.  Its hand-written gradient includes the derivatives of `S_full`,
`B1`, `Z`, `P`, every diagonal shell norm, and every normalized absolute Gram
entry in the chosen row.  Running all rows separately covers the nonsmooth
finite maximum `R_H`.  Initial/final central-difference errors in the current
K8 certificates are at most `1.66e-10`; the final K12 errors are at most
`1.50e-10`.

```text
branch                         initial K8    K8 16-step   K12 zero-pad  K12 1-step
H8, row shell 5              3.44363e-4    3.63935e-4   3.68341e-4    3.68688e-4
H16, row shell 6             3.71126e-4    3.97807e-4   4.12917e-4    4.18215e-4
```

The H16 K12 result is the current finite stress record.  It is `7.46%` above
the previous `3.89149e-4` maximum inherited from Cauchy-majorant optimization.
On the H16 K8 branch, `D_H` rises by `9.30%` and the chosen Gram row by
`2.56%`, while the common normalization falls by `4.37%`; the joint objective
still rises by `7.19%`.  The K12 step repeats the same competition and gains
another `1.283%`.  This is direct evidence that optimizing a proxy left room
in PNT-12.  It neither proves growth without bound nor supplies a uniform
constant.  An independent full Gram ledger confirms that shell 6 remains the
maximal row (`R_H=1.52102`) and reconstructs the final `4.18215e-4` value
exactly.

To keep this exact search within workstation memory, the triad engine now has
a direct `vjp_sum` kernel over group-index lists.  It avoids materializing
copied aggregate interaction families.  On the same K12 check, peak RSS fell
from `7.94 GiB` to `4.88 GiB` with the objective unchanged; deterministic
`static,1` shell scheduling also makes the saved L-BFGS trace replayable.

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
  --line-search 14 --lbfgs-history 8 --threads 12 --step 0.1

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

./build/navier_stokes_lab local-sld-projective-normalization-cauchy-scan \
  --height 8 --state proof/l4/states/local-projective-normalization-cauchy/H8-K12-one-step/K12.tsv \
  --height 16 --state proof/l4/states/local-projective-normalization-cauchy/H16-K12-one-step/K12.tsv \
  --height 32 --state proof/l4/states/local-projective-normalization-cauchy/H32-K12-one-step/K12.tsv \
  --height 64 --state proof/l4/states/local-projective-normalization-cauchy/H64-K12-one-step/K12.tsv \
  --selection double-triple-remainder-without-123 --threads 12 \
  --certificate /tmp/K12-normalization-cauchy-height-scan.json

./build/navier_stokes_lab local-sld-projective-normalization-satellite \
  --base-state proof/l4/states/local-projective-normalization-cauchy/H64-K8-twenty-four-step/K8.tsv \
  --base-cutoff 2 --cutoff 5 --cutoff 6 --cutoff 8 --cutoff 10 --cutoff 12 \
  --projective-core-height 64 --satellite-coefficient 1 \
  --satellite-modes 0 --selection double-triple-remainder-without-123 \
  --seed 20260802 --threads 12 --state-dir /tmp/normalization-satellite \
  --certificate /tmp/normalization-satellite.json

./build/navier_stokes_lab local-sld-projective-normalization-tail-schur \
  --height 8 --state proof/l4/states/local-projective-normalization-cauchy/H8-K12-one-step/K12.tsv \
  --height 16 --state proof/l4/states/local-projective-normalization-cauchy/H16-K12-one-step/K12.tsv \
  --height 32 --state proof/l4/states/local-projective-normalization-cauchy/H32-K12-one-step/K12.tsv \
  --height 64 --state proof/l4/states/local-projective-normalization-cauchy/H64-K12-one-step/K12.tsv \
  --selection double-triple-remainder-without-123 --threads 12 \
  --certificate /tmp/K12-normalization-tail-schur.json

./build/navier_stokes_lab local-sld-projective-height-gap-correlation \
  --state proof/l4/states/local-projective-normalization-cauchy/H8-K8-twenty-four-step/K8.tsv \
  --output-state /tmp/K8-gap6.tsv --certificate /tmp/K8-gap6.json \
  --selection double-triple-remainder-without-123 \
  --first-shell 2 --second-shell 8 --iterations 24 \
  --line-search 14 --lbfgs-history 8 --threads 12 --step 0.1

./build/navier_stokes_lab local-sld-projective-height-gap-output \
  --state proof/l4/states/local-projective-height-gap-correlation/K8-shell2-shell8-twenty-four-step/K8.tsv \
  --selection double-triple-remainder-without-123 \
  --first-shell 2 --second-shell 8 --top-modes 32 --threads 12 \
  --certificate /tmp/K8-gap6-output.json

./build/navier_stokes_lab local-sld-projective-height-gap-triad-attribution \
  --state proof/l4/states/local-projective-height-gap-correlation/K8-shell2-shell8-twenty-four-step/K8.tsv \
  --selection double-triple-remainder-without-123 \
  --first-shell 2 --second-shell 8 \
  --output-x 8 --output-y -1 --output-z 8 \
  --top-interactions 32 --threads 12 \
  --certificate /tmp/K8-gap6-triads.json

./build/navier_stokes_lab local-sld-projective-normalization-schur \
  --state proof/l4/states/local-projective-normalization-cauchy/H16-K8-twenty-four-step/K8.tsv \
  --output-state /tmp/H16-row6-K8.tsv \
  --certificate /tmp/H16-row6-K8.json \
  --selection double-triple-remainder-without-123 \
  --projective-core-height 16 --row-shell 6 \
  --iterations 16 --line-search 14 --lbfgs-history 8 \
  --threads 12 --step 0.1
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
- [`../../analysis/shifted-local-density/remainder-quartet/K8-normalization-cauchy-H8-H64-lightweight-scan.json`](../../analysis/shifted-local-density/remainder-quartet/K8-normalization-cauchy-H8-H64-lightweight-scan.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K12-zero-pad-normalization-cauchy-H8-H64-lightweight-scan.json`](../../analysis/shifted-local-density/remainder-quartet/K12-zero-pad-normalization-cauchy-H8-H64-lightweight-scan.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K12-one-step-normalization-cauchy-H8-H64-lightweight-scan.json`](../../analysis/shifted-local-density/remainder-quartet/K12-one-step-normalization-cauchy-H8-H64-lightweight-scan.json)
- [`../../adversary/shifted-local-density/projective-normalization-cauchy/H128-K12-four-step.json`](../../adversary/shifted-local-density/projective-normalization-cauchy/H128-K12-four-step.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K4-K12-normalization-satellite-H64-baseK2-c1-m1.json`](../../analysis/shifted-local-density/remainder-quartet/K4-K12-normalization-satellite-H64-baseK2-c1-m1.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K5-K12-normalization-satellite-H64-baseK2-c1-dense.json`](../../analysis/shifted-local-density/remainder-quartet/K5-K12-normalization-satellite-H64-baseK2-c1-dense.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K8-normalization-cauchy-H8-H64-tail-schur.json`](../../analysis/shifted-local-density/remainder-quartet/K8-normalization-cauchy-H8-H64-tail-schur.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K12-one-step-normalization-cauchy-H8-H64-tail-schur.json`](../../analysis/shifted-local-density/remainder-quartet/K12-one-step-normalization-cauchy-H8-H64-tail-schur.json)
- [`../../adversary/shifted-local-density/projective-height-gap-correlation/K8-shell2-shell8-twenty-four-step.json`](../../adversary/shifted-local-density/projective-height-gap-correlation/K8-shell2-shell8-twenty-four-step.json)
- [`../../adversary/shifted-local-density/projective-height-gap-correlation/K12-shell2-shell9-twenty-four-step.json`](../../adversary/shifted-local-density/projective-height-gap-correlation/K12-shell2-shell9-twenty-four-step.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K8-K12-height-gap-correlation-PNT12-tradeoff-tail-schur.json`](../../analysis/shifted-local-density/remainder-quartet/K8-K12-height-gap-correlation-PNT12-tradeoff-tail-schur.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K8-shell2-shell8-height-gap-shared-output.json`](../../analysis/shifted-local-density/remainder-quartet/K8-shell2-shell8-height-gap-shared-output.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K12-shell2-shell9-height-gap-shared-output.json`](../../analysis/shifted-local-density/remainder-quartet/K12-shell2-shell9-height-gap-shared-output.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K8-shell2-shell8-dominant-output-triad-attribution.json`](../../analysis/shifted-local-density/remainder-quartet/K8-shell2-shell8-dominant-output-triad-attribution.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K12-shell2-shell9-dominant-output-triad-attribution.json`](../../analysis/shifted-local-density/remainder-quartet/K12-shell2-shell9-dominant-output-triad-attribution.json)
- [`../../adversary/shifted-local-density/projective-normalization-schur/H8-row5-K8-sixteen-step.json`](../../adversary/shifted-local-density/projective-normalization-schur/H8-row5-K8-sixteen-step.json)
- [`../../adversary/shifted-local-density/projective-normalization-schur/H8-row5-K12-one-step.json`](../../adversary/shifted-local-density/projective-normalization-schur/H8-row5-K12-one-step.json)
- [`../../adversary/shifted-local-density/projective-normalization-schur/H16-row6-K8-sixteen-step.json`](../../adversary/shifted-local-density/projective-normalization-schur/H16-row6-K8-sixteen-step.json)
- [`../../analysis/shifted-local-density/remainder-quartet/H16-row6-K8-sixteen-step-zero-pad-K12-tail-schur.json`](../../analysis/shifted-local-density/remainder-quartet/H16-row6-K8-sixteen-step-zero-pad-K12-tail-schur.json)
- [`../../adversary/shifted-local-density/projective-normalization-schur/H16-row6-K12-one-step.json`](../../adversary/shifted-local-density/projective-normalization-schur/H16-row6-K12-one-step.json)
- [`../../analysis/shifted-local-density/remainder-quartet/H16-row6-K12-one-step-tail-schur.json`](../../analysis/shifted-local-density/remainder-quartet/H16-row6-K12-one-step-tail-schur.json)
