# Dyadic projective-height Schur reduction

This note gives an exact algebraic reduction of the open projective part of
RQ-11 to one cutoff-uniform joint dyadic estimate. It does not prove that
uniform estimate and therefore does not prove RQ-11, the local SLD lemma, or
the Clay problem.

## Exact height matrix

Partition the primitive feasible squared-length shapes by

```text
Gamma_0={c=1},
Gamma_j={2^(j-1)<c<=2^j},  j>=1,
```

where `c` is the largest primitive squared length. Write `B_j` for the sum of
all projective advection operators in `Gamma_j`, and set

```text
b_j=B_j(u,u),  c_j=B_j(u,Au),
S_j=<Au,b_j>,  T_j=<A b_j,Au>.
```

The diagonal quartet `J_jj` and the symmetric off-diagonal quartet `J_ij`
are the five-component polarization of the closed bracket:

```text
J_jj = -<b_j,A b_j> + <b_j,c_j> - <Au,B_j(b_j,u)>
       + S_j^2/(2Z) + 3 S_j T_j/(2P),

J_ij = -<b_i,A b_j>-<b_j,A b_i>
       +<b_i,c_j>+<b_j,c_i>
       -<Au,B_i(b_j,u)>-<Au,B_j(b_i,u)>
       +S_i S_j/Z+3(S_i T_j+S_j T_i)/(2P).          (DHS-1)
```

`LocalSldProjectiveHeightMatrix` evaluates every term in DHS-1 directly. Its
sum reconstructs the selected full bracket below `1e-17` relative error on
all saved stress states.

## Absolute-component Schur lemma

Let

```text
q=|S_full|/(Z^2 P^2).
```

Define `e_j` as `q` times the sum of the absolute values of the five
components of `J_jj`, and define `e_ij` analogously from the five components
of `J_ij`. For every pair with nonzero diagonal weights set

```text
r_jj=1,
r_ij=e_ij/(2 sqrt(e_i e_j)),  i!=j,
R=max_i sum_j r_ij.                                  (DHS-2)
```

Then the arithmetic--geometric mean inequality gives the exact finite-matrix
bound

```text
sum_j e_j + sum_(i<j) e_ij
 <= sum_j e_j + sum_(i<j) r_ij(e_i+e_j)
 <= R sum_j e_j.                                    (DHS-3)
```

The left side dominates the absolute value of the complete projective
power-one block. Therefore the direct sufficient condition for RQ-11 is

```text
sup_(N,u) R_N(u) sum_j e_j,N(u) < infinity.          (DHS-4)
```

The two separate estimates

```text
sup_(N,u) R_N(u) < infinity,
sup_(N,u) sum_j e_j,N(u) < infinity                  (DHS-5)
```

would imply DHS-4, but they are stronger than necessary. A sufficient form
of the second estimate is a dyadic decay `e_j<=C 2^(-alpha j)` for any
`alpha>0`. DHS-3 is proved algebraically; the joint estimate DHS-4 is the
remaining analytic statement.

The component envelope is essential. Normalizing off-diagonal entries by the
absolute value of the signed diagonal `J_jj` is invalid because internal
cancellation can make `J_jj` almost zero. The finite ratio then exceeds
`10^7` on a saved K8 state. The five-component envelope removes that
artificial small denominator.

## Current falsification results

On the H=64 K8 open-power winner, the exact matrix gives

```text
total component envelope      0.0038520661
diagonal component envelope   0.0021599994
maximum Schur row sum          2.9859137
Schur upper bound              0.0064495718
actual / Schur upper bound     0.5972592.
```

The maximum normalized off-diagonal envelope is `0.5372`. By shell gap its
observed maxima decrease as

```text
gap       1       2       3       4       5       6       7       8
ratio   .537    .426    .283    .246    .158    .0529   .00789  .000125.
```

The maximum row sum lies between `2.89` and `3.85` on the initial K8
open-power, height-power, and height-stretching adversary winners. A coherent
fan perturbed by the same K8 state gives row sums `2.987`, `3.993`, and
`4.089` at cutoffs K8, K12, and K16. Its corresponding Schur products
`R sum_j e_j` are `1.55e-10`, `5.45e-12`, and `3.79e-13`. Thus a separately
uniform row-sum estimate may be unnecessarily strong: the row factor grows
on this finite sequence while the joint quantity in DHS-4 collapses. This is
finite diagnostic evidence, not a theorem that either separate bound fails.

An exact-gradient objective maximizes the complete signed diagonal block

```text
|J_(H,2H] S_full|^2/(Z^4P^4).
```

On separately optimized K8 states its square roots at H=8,16,32,64 are
approximately

```text
5.47e-4, 4.73e-4, 3.93e-4, 2.45e-4.
```

The fitted signed-diagonal height slope is `-0.374`. More importantly, the
absolute five-component diagonal envelope has fitted slope `-0.435`, so the
observed decay is not caused only by cancellation inside `J_jj`. These are
optimized finite lower branches, not an upper bound in DHS-4 or DHS-5.

The exact-gradient `projective-height-outer-power-ratio` objective attacks the
leading coercive part of the joint quantity directly:

```text
|S_full| sum_j ||A^(1/2) B_j(u,u)||_2^2/(Z^2 P^2).
```

Starting from the H=64 K8 height-power winner, 36 L-BFGS iterations in two
restarted stages raise this quantity from `0.001472` to `0.002940`. The exact
full component matrix of the resulting state has `R=2.8952`, diagonal
envelope `0.0051801`, total envelope `0.0089833`, and Schur upper bound
`0.0149976`.

The sharper exact-gradient `projective-height-envelope-ratio` objective now
maximizes the left side of DHS-3 itself:

```text
[ |S_full| sum_(i<=j) sum_(five components) |J_ij,component|
  /(Z^2 P^2) ]^2.                                  (DHS-6)
```

Its sign-chamber reverse derivative includes all shell pairs and all five
components; a central-difference test has relative error `6.61e-12`. Eight
K8 L-BFGS steps raise the unsquared envelope from `0.0089833` to `0.0094401`.
The resulting exact matrix has diagonal envelope `0.0052205`, `R=2.9616`,
and Schur bound `0.0154608`.

A sparse-support K12 continuation reaches `0.0094648` after six steps. Sparse
support omits inactive Fourier targets and is explicitly marked
`complete_galerkin_cutoff=false`; it is only a restricted-support diagnostic.
Zero-padding that winner to the complete K12 Galerkin cube gives `0.0095014`
without any further K12 optimization. Its exact matrix has diagonal envelope
`0.0052331`, `R=2.9681`, and Schur bound `0.0155326`. Thus the optimized finite
total envelope does not exhibit the earlier diagonal-shell decay between K8
and K12. This neither disproves nor proves a cutoff-uniform bound; it removes
finite decay as the current reason to expect one.

## Reproduction

```bash
./build/navier_stokes_lab local-sld-projective-height-matrix \
  --state proof/l4/states/local-projective-open-power/height-64-K8-warm/K8.tsv \
  --threads 12 --exclude-triple-family \
  --certificate proof/l4/analysis/shifted-local-density/remainder-quartet/K8-open-power-height64-winner-height-matrix.json

./build/navier_stokes_lab local-sld-projective-height-transfer \
  --height 8 --state proof/l4/states/local-projective-height-power/H8-K8-warm/K8.tsv \
  --height 16 --state proof/l4/states/local-projective-height-power/H16-K8-warm/K8.tsv \
  --height 32 --state proof/l4/states/local-projective-height-power/H32-K8-warm/K8.tsv \
  --height 64 --state proof/l4/states/local-projective-height-power/H64-K8-warm/K8.tsv \
  --threads 12 --exclude-triple-family \
  --certificate proof/l4/analysis/shifted-local-density/remainder-quartet/K8-optimized-height-power-transfer.json

./build/navier_stokes_lab local-closure-adversary \
  --objective projective-height-outer-power-ratio \
  --selection double-triple-remainder --min-cutoff 8 --max-cutoff 8 \
  --restarts 1 --workers 12 --iterations 12 --method lbfgs --backend direct \
  --warm-state proof/l4/states/local-projective-height-power/H64-K8-warm/K8.tsv \
  --certificate proof/l4/adversary/shifted-local-density/projective-height-outer-power/K8-height-power-warm.json \
  --state-dir proof/l4/states/local-projective-height-outer-power/K8-height-power-warm

./build/navier_stokes_lab local-closure-adversary \
  --objective projective-height-envelope-ratio \
  --selection double-triple-remainder --min-cutoff 12 --max-cutoff 12 \
  --restarts 1 --workers 12 --iterations 4 --method lbfgs --backend direct \
  --lean --preserve-warm-layout \
  --warm-state proof/l4/states/local-projective-height-envelope/K12-warm-blend-two-step-sparse/K12.tsv \
  --certificate proof/l4/adversary/shifted-local-density/projective-height-envelope/K12-six-step-sparse.json \
  --state-dir proof/l4/states/local-projective-height-envelope/K12-six-step-sparse
```

Artifacts:

- [`../../analysis/shifted-local-density/remainder-quartet/K8-open-power-height64-winner-height-matrix.json`](../../analysis/shifted-local-density/remainder-quartet/K8-open-power-height64-winner-height-matrix.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K8-optimized-height-power-transfer.json`](../../analysis/shifted-local-density/remainder-quartet/K8-optimized-height-power-transfer.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K8-height-outer-power-winner-height-matrix.json`](../../analysis/shifted-local-density/remainder-quartet/K8-height-outer-power-winner-height-matrix.json)
- [`../../adversary/shifted-local-density/projective-height-outer-power/K8-height-power-warm.json`](../../adversary/shifted-local-density/projective-height-outer-power/K8-height-power-warm.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K8-height-envelope-winner-height-matrix.json`](../../analysis/shifted-local-density/remainder-quartet/K8-height-envelope-winner-height-matrix.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K12-height-envelope-six-step-full-matrix.json`](../../analysis/shifted-local-density/remainder-quartet/K12-height-envelope-six-step-full-matrix.json)
- [`../../adversary/shifted-local-density/projective-height-envelope/K12-six-step-sparse.json`](../../adversary/shifted-local-density/projective-height-envelope/K12-six-step-sparse.json)
- [`../../adversary/shifted-local-density/projective-height-envelope/K12-six-step-full-evaluate.json`](../../adversary/shifted-local-density/projective-height-envelope/K12-six-step-full-evaluate.json)
