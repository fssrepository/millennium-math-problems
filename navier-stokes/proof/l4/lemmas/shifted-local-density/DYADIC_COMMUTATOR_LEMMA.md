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

DCL-6 and DCL-7 remain open. They are stronger than the joint condition
DCL-5, but they separate the derivative-bearing commutator geometry from the
outer coercive weight.

## Why the weight must be coercive

Using the paired diagonal `g_jj` itself as the Schur weight fails even on a
finite stress state. The height-one diagonal is about `1e-25` because DCL-1
cancels, while its off-diagonal couplings are nonzero. The resulting row sum
exceeds `9e7`. This is a normalization failure, not a large quartet.

The outer weight `w_j` does not contain that cancellation. It gives no
unscaled pair on the saved K8 and K12 complete-Galerkin winners.

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
