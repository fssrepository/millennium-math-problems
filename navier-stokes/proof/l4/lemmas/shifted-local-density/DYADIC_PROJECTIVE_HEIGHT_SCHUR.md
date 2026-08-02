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

## Commutator-paired coercive Schur lemma

The five-component envelope is sufficient but unnecessarily separates its
two largest terms. Preserve their exact algebraic cancellation by defining

```text
g_ij=q( |outer_ij+advected_ij|+|nested_ij|
          +|enstrophy_ij|+|palinstrophy_ij| ).       (DHS-6)
```

The signed quartet is still bounded by `sum_(i<=j) g_ij`. Define the
coercive outer weights

```text
w_i=q |outer_ii|=q ||A^(1/2)B_i(u,u)||_2^2.         (DHS-7)
```

Unlike the signed or commutator-paired diagonal, `w_i` does not lose scale
through internal cancellation. Set

```text
rho_ii=g_ii/w_i,
rho_ij=g_ij/(2 sqrt(w_i w_j)), i!=j,
R_comm=max_i sum_j rho_ij.                           (DHS-8)
```

The same arithmetic--geometric mean argument proves the exact finite bound

```text
sum_(i<=j) g_ij <= R_comm sum_i w_i.                 (DHS-9)
```

Consequently the current sharp sufficient condition for RQ-11 is

```text
sup_(N,u) R_comm,N(u) sum_i w_i,N(u) < infinity.     (DHS-10)
```

A natural analytic route is an off-diagonal estimate for `rho_ij` that
decays in `|i-j|`, together with a bound for the outer weight. Neither is
proved here. Normalizing `g_ij` by `sqrt(g_ii g_jj)` is invalid: on a saved
K8 winner shell zero has paired diagonal about `1e-25` but nonzero coupling
to higher shells, producing a spurious row sum above `9e7`. The outer weights
in DHS-7 remove exactly this cancellation defect.

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
  /(Z^2 P^2) ]^2.                                  (DHS-11)
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

The commutator-paired majorant in DHS-6 has its own exact-gradient objective,
`projective-height-commutator-envelope-ratio`, with central-difference error
`7.74e-12`. Direct optimization gives the following complete-Galerkin stress
data:

```text
                         K8              K12
paired envelope       0.00179960       0.00182423
signed block          0.00082522       0.00081183
outer weight          0.00273347       0.00275855
R_comm                0.881010         0.883351
Schur upper bound     0.00240822       0.00243677.
```

The finite inequality DHS-9 is verified in both certificates and has no
unscaled off-diagonal pair. The paired envelope is about five times below
the separately absolute five-component envelope and only about twice the
signed target. Its slight K8-to-K12 increase rules out finite decay as
evidence, but is fully compatible with a cutoff-uniform constant.

The separate row-bound route is not the active target anymore. An exact K8
gradient search raises the commutator-paired envelope divided by the outer
weight from `0.6584` to `4.8496` while moving toward a small-outer-weight
state. Pairing the nested term as well gives the exact dynamic entry
`D_ij=outer_ij+advected_ij+nested_ij`; its isolated outer-weight ratio is also
ill-conditioned. The relevant joint objective, however, remains finite in
the current search: the dynamic-paired normalized envelope reaches
`0.00186662` after 44 K8 steps, and its complete zero-padded K12 value is
`0.00187172`; one full K12 step raises it to `0.00187877`. The K8 exact matrix verifies its finite outer Schur inequality
with row sum `0.857683`, bound `0.00248948`, and actual/bound ratio `0.749804`.
The corresponding one-step K12 values are `0.875077`, `0.00254378`, and
`0.738575`. Its maximum finite `2^gap rho_gap` diagnostic is `0.63863`.
None of these finite values proves the required cutoff-uniform joint estimate.

The dynamic entry also has an exact response representation. With
`R_j=C_j-B_j(.,u)^*Au`, it is `<b_j,R_j>` on the diagonal and the symmetric
cross pairing off the diagonal. Replacing the outer-only weight by
`q(||A^(1/2)b_j||^2+||A^(-1/2)R_j||^2)` reduces the finite Schur row from
`0.857683` to `0.672016` at K8 and from `0.875077` to `0.684904` at K12. On
the outer-degenerate stress state it reduces `4.55082` to `1.14162`. This
removes the observed normalization singularity but does not prove a uniform
response-weighted Schur product.

For the signed dynamical block required by RQ-11, a stronger simplification is
available. Summing the exact response entries before taking absolute values
gives `<sum_j b_j,sum_j R_j>`. Its absolute value is bounded directly by the
global `H1`--`H-1` Young weight, so no cutoff-uniform Schur row or gap-decay
theorem is needed for that signed block. The K8 reconstruction error is
`1.37e-18`; its signed normalized dynamical value `0.00136712` is below the
direct global bound `0.00298882`. The height-Schur data remain useful stress
diagnostics and for stronger shellwise envelopes, but the active proof target
has narrowed to bounding the global response weight and the two exact
normalization terms.

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

./build/navier_stokes_lab local-closure-adversary \
  --objective projective-height-commutator-envelope-ratio \
  --selection double-triple-remainder --min-cutoff 8 --max-cutoff 8 \
  --restarts 1 --workers 12 --iterations 12 --method lbfgs --backend direct \
  --lean \
  --warm-state proof/l4/states/local-projective-height-commutator-envelope/K8-eight-step/K8.tsv \
  --certificate proof/l4/adversary/shifted-local-density/projective-height-commutator-envelope/K8-twenty-step.json \
  --state-dir proof/l4/states/local-projective-height-commutator-envelope/K8-twenty-step
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
- [`../../analysis/shifted-local-density/remainder-quartet/K8-height-commutator-envelope-twenty-step-matrix.json`](../../analysis/shifted-local-density/remainder-quartet/K8-height-commutator-envelope-twenty-step-matrix.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K12-height-commutator-envelope-twelve-step-full-matrix.json`](../../analysis/shifted-local-density/remainder-quartet/K12-height-commutator-envelope-twelve-step-full-matrix.json)
- [`../../adversary/shifted-local-density/projective-height-commutator-envelope/K8-twenty-step.json`](../../adversary/shifted-local-density/projective-height-commutator-envelope/K8-twenty-step.json)
- [`../../adversary/shifted-local-density/projective-height-commutator-envelope/K12-twelve-step-full-evaluate.json`](../../adversary/shifted-local-density/projective-height-commutator-envelope/K12-twelve-step-full-evaluate.json)
