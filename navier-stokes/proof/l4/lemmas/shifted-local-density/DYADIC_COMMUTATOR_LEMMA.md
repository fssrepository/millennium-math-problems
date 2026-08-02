# Dyadic projective commutator target

This note isolates the current lowest analytic target inside the open
projective part of RQ-11. It proves the algebraic reduction, not the required
cutoff-uniform estimate. It therefore does not prove the local SLD lemma or
the Clay problem.

## Exact commutator identity

Let `B_j` be the projective advection operator in the `j`-th dyadic primitive
height shell and write

```text
b_j=B_j(u,u),
t_j=B_j(u,Au),
C_j=t_j-A b_j.
```

The sum of the outer-square and advected components of the height quartet is

```text
outer_jj+advected_jj = <b_j,C_j>,                   (DCL-1)

outer_ij+advected_ij = <b_i,C_j>+<b_j,C_i>, i!=j.  (DCL-2)
```

These identities follow directly from self-adjointness of `A`. In Fourier
variables the commutator is

```text
C_j(k)=P_k sum_(p+q=k, shape(p,q,k) in Gamma_j)
          i(q dot u_p)(|q|^2-|k|^2)u_q.             (DCL-3)
```

Thus the symbol difference that carries the derivative gain is retained
before any absolute value is taken.

## Coercive Schur reduction

Set

```text
q=|S_full|/(Z^2P^2),
w_j=q ||A^(1/2)b_j||_2^2,

g_ij=q( |outer_ij+advected_ij|+|nested_ij|
          +|enstrophy_ij|+|palinstrophy_ij| ).
```

For nonzero weights define

```text
rho_jj=g_jj/w_j,
rho_ij=g_ij/(2 sqrt(w_iw_j)), i!=j,
R_comm=max_i sum_j rho_ij.
```

Arithmetic--geometric mean gives the exact finite inequality

```text
sum_(i<=j) g_ij <= R_comm sum_j w_j.                (DCL-4)
```

The left side dominates the signed projective power-one block. Hence RQ-11
would follow from

```text
sup_(N,u) R_comm,N(u) sum_j w_j,N(u) < infinity.    (DCL-5)
```

Two sufficient analytic estimates are

```text
rho_ij <= C 2^(-alpha |i-j|), alpha>0,              (DCL-6)
sup_(N,u) sum_j w_j,N(u) < infinity.                (DCL-7)
```

DCL-6 and DCL-7 are stronger than the joint condition DCL-5. Exact-gradient
stress tests now show that this separation is badly conditioned: at fixed K8,
maximizing `(sum g_ij)/(sum w_j)` raises the ratio from `0.6584` to `4.8496`
while both numerator and denominator decrease. The projected gradient remains
large. This is finite numerical evidence, not an analytic counterexample, but
it removes a cutoff-uniform bound on `R_comm` as the active proof target.

## Full dynamic pairing

The same stress branch identifies a second exact cancellation. Define the
complete dynamical entry

```text
D_ij=outer_ij+advected_ij+nested_ij.                (DCL-8)
```

On the strongest outer-coercivity stress state, the globally signed
commutator and nested contributions are `+4.83e-8` and `-1.07e-7` after
power-one normalization. Taking their absolute values separately discards
this cancellation. The sharper matrix majorant is therefore

```text
h_ij=q( |D_ij|+|enstrophy_ij|+|palinstrophy_ij| ). (DCL-9)
```

The code computes DCL-9 directly and differentiates its sign chambers
exactly. With the same outer weights, finite-dimensional AM--GM gives the
corresponding exact Schur bound. The ratio `(sum h_ij)/(sum w_j)` can also be
driven upward at fixed K8 (`1.823` to `3.209` in 24 steps), so separating a
uniform row constant from the vanishing outer weight remains too strong.

The active joint target is instead the scale-normalized quantity itself:

```text
sup_(N,u) sum_(i<=j) h_ij,N(u) < infinity.          (DCL-10)
```

Unlike the isolated row constant, DCL-10 is exactly the type of joint bound
needed by RQ-11. Forty-four K8 exact-gradient steps raise its unsquared value
from `0.00177488` to `0.00186662`. Zero-padding that state to the complete K12
Galerkin cube gives `0.00187172`; one full K12 gradient step raises it to
`0.00187877`. The activated gap-nine contribution is `3.92e-16`, and the
finite diagnostic satisfies `rho_gap <= 0.63863 2^(-gap)`. These are finite
lower branches and fitted diagnostics, not an upper bound or a proof of
DCL-10.

## Dynamic response weight

The outer-only normalization fails near states for which `b_j` is small but
the linear response responsible for the nested term is not comparably small.
Let `B_j(.,u)^*` denote the adjoint in the advecting slot and define

```text
N_j=-B_j(.,u)^* Au,
R_j=C_j+N_j.                                        (DCL-11)
```

The full dynamical matrix then has the exact mixed-Gram representation

```text
D_jj=<b_j,R_j>,
D_ij=<b_i,R_j>+<b_j,R_i>, i!=j.                    (DCL-12)
```

The C++ matrix engine evaluates `R_j` with the exact bilinear VJP and checks
DCL-12 entry by entry. The maximum reconstruction errors are `2.17e-16` at
K8 and `1.44e-15` at K12. This suggests the nondegenerate combined weights

```text
v_j=q( ||A^(1/2)b_j||_2^2+||A^(-1/2)R_j||_2^2 ).  (DCL-13)
```

Normalizing DCL-9 by `v_j` gives another exact finite Schur inequality. On the
joint-envelope winners its `(weight,row,bound)` values are

```text
K8   (0.00359195, 0.672016, 0.00241385)
K12  (0.00359895, 0.684904, 0.00246494).
```

On the outer-degenerate K8 stress state, the outer row is `4.55082`, whereas
the response-weighted row is `1.14162`; no pair is unscaled. Thus DCL-13 fixes
the specific degeneration found by the adversary. The remaining proof task is
still cutoff-uniform: bound the response-weighted Schur product analytically.
The finite data do not establish such a bound.

## Global response reduction

The shellwise envelope is stronger than the signed statement needed by
RQ-11. Set `b=sum_j b_j` and `R=sum_j R_j`. Summing DCL-12 first gives the
exact identity

```text
sum_(i<=j) D_ij=<b,R>.                             (DCL-14)
```

Therefore Young's inequality in the `H1`--`H-1` duality gives, without a
Schur row or a shell-gap assumption,

```text
q |sum_(i<=j) D_ij|
 <= q( ||A^(1/2)b||_2^2+||A^(-1/2)R||_2^2 )/2.   (DCL-15)
```

The matrix engine reconstructs DCL-14 from all entries. On the K8 joint
winner its relative error is `1.37e-18`; the signed dynamical value is
`0.00136712`, while the DCL-15 bound is `0.00298882` with ratio `0.457413`.
On the outer-degenerate stress state the corresponding bound ratio is
`0.619884`. Thus Schur summability is no longer required for the signed
dynamical part of RQ-11.

Writing `s=<Au,b>` and `t=<Ab,Au>`, the complete selected bracket is exactly

```text
<b,R> + s^2/(2Z) + 3st/(2P).                       (DCL-16)
```

The first normalization is bounded by `||A^(1/2)b||^2/2`. Direct Cauchy gives
`3 sqrt(Z/P) ||A^(1/2)b|| ||Ab||/2` for the second, but this bound is loose on
the current winner. The remaining analytic task is a cutoff-uniform estimate
for the DCL-15 response weight together with these two global normalization
terms, especially the palinstrophy cross term. No finite certificate proves
that estimate.

## Why the paired diagonal is not a weight

Using the paired diagonal `g_jj` itself as the Schur weight fails even on a
finite stress state. The height-one diagonal is about `1e-25` because DCL-1
cancels, while its off-diagonal couplings are nonzero. The resulting row sum
exceeds `9e7`. This is a normalization failure, not a large quartet.

The outer weight `w_j` does not contain that diagonal cancellation and gives
no unscaled pair on the saved K8 and K12 complete-Galerkin winners. DCL-13
adds the missing dynamic response needed near the outer null set.

## Current exact diagnostics

```text
                         K8              K12
paired envelope       0.00179960       0.00182423
signed block          0.00082522       0.00081183
sum_j w_j             0.00273347       0.00275855
R_comm                0.881010         0.883351
DCL-4 upper bound     0.00240822       0.00243677
actual / bound        0.747275         0.748627.
```

The maximum observed `rho_ij` by shell gap is

```text
gap       0       1       2       3       4       5        6         7          8
K8     .6245   .1314   .0908   .0753   .0387   .0176   .00186   .000693   .0000156
K12    .6622   .1315   .0916   .0759   .0393   .0175   .00185   .000709   .00000554
```

K12 also has gap-nine maximum `4.69e-7`. This is finite evidence for the
form of DCL-6, not an upper bound or a fitted proof.

Artifacts:

- [`../../analysis/shifted-local-density/remainder-quartet/K8-height-commutator-envelope-twenty-step-matrix.json`](../../analysis/shifted-local-density/remainder-quartet/K8-height-commutator-envelope-twenty-step-matrix.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K12-height-commutator-envelope-twelve-step-full-matrix.json`](../../analysis/shifted-local-density/remainder-quartet/K12-height-commutator-envelope-twelve-step-full-matrix.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K8-height-dynamic-envelope-twenty-four-step-matrix.json`](../../analysis/shifted-local-density/remainder-quartet/K8-height-dynamic-envelope-twenty-four-step-matrix.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K8-height-dynamic-envelope-forty-four-step-matrix.json`](../../analysis/shifted-local-density/remainder-quartet/K8-height-dynamic-envelope-forty-four-step-matrix.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K12-height-dynamic-envelope-forty-four-step-zero-pad-matrix.json`](../../analysis/shifted-local-density/remainder-quartet/K12-height-dynamic-envelope-forty-four-step-zero-pad-matrix.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K12-height-dynamic-envelope-one-step-full-matrix.json`](../../analysis/shifted-local-density/remainder-quartet/K12-height-dynamic-envelope-one-step-full-matrix.json)
- [`../../adversary/shifted-local-density/projective-height-dynamic-envelope/K8-twenty-four-step.json`](../../adversary/shifted-local-density/projective-height-dynamic-envelope/K8-twenty-four-step.json)
- [`../../adversary/shifted-local-density/projective-height-dynamic-envelope/K8-forty-four-step.json`](../../adversary/shifted-local-density/projective-height-dynamic-envelope/K8-forty-four-step.json)
- [`../../adversary/shifted-local-density/projective-height-dynamic-envelope/K12-forty-four-step-zero-pad-evaluate.json`](../../adversary/shifted-local-density/projective-height-dynamic-envelope/K12-forty-four-step-zero-pad-evaluate.json)
- [`../../adversary/shifted-local-density/projective-height-dynamic-envelope/K12-one-step-full.json`](../../adversary/shifted-local-density/projective-height-dynamic-envelope/K12-one-step-full.json)
