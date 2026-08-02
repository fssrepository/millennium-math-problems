# Weighted response interaction tensor

This note records an explicit cutoff-diagonal bilinear candidate extracted
from the exact direct-triad interaction tensor. It is a finite computational
certificate and does not prove the candidate, the shifted local-density lemma,
or the Clay problem.

## Exact finite tensor

Let `b_0,...,b_K` be the orthonormal boundary-free cyclic response basis from
`RESPONSE_DIAGONAL.md`. The code evaluates every ordered interaction

```text
T_mij = <b_m, B(b_i,b_j)>,              0 <= i,j,m <= K,
```

with direct Fourier-triad summation. For input radius `R` and output radius
`r`, it reports

```text
C_K(R,r) = max_(i,j) [sum_(m=0)^K r^m |T_mij|] / R^(i+j).   (RT-1)
```

Using separate radii is essential. With `R=r=1.25`, `C_K` grows from
`0.9723` at K2 to `5.2028` at K8. This is the expected one-derivative loss;
the same-radius weighted algebra is not a viable closure.

## Explicit surviving candidate

For `R=2` and `r=23/20=1.15`, all ordered pairs through K12 satisfy

```text
C_K(2,23/20) = 23/(20 sqrt(3))
              = 0.6639528095680697,                         (RT-2)
```

and equality is attained by the axis interaction `(i,j)=(0,0)`. This gives
the concrete infinite-dimensional candidate

```text
RT-L:
sum_(m>=0) (23/20)^m |<b_m,B(b_i,b_j)>|
    <= [23/(20 sqrt(3))] 2^(i+j)       for every i,j >= 0.  (RT-L)
```

The equality constant is not fitted: the exact first interaction has
`|<b_1,B(b_0,b_0)>|=1/sqrt(3)`. The machine result says that no later
boundary-free pair through K12 exceeds this analytic candidate. It does not
establish RT-L for arbitrary response order.

At the less separated output radius `r=1.25`, the measured maxima remain
bounded by `0.962525` through K12 but no longer stay at the axis interaction.
The wider `R=2 -> r=1.15` gap exposes the simpler sharp structure.

## Orthogonal complement

For every pair the tensor engine also computes

```text
w_perp(i,j) = ||(I-P_K)B(b_i,b_j)||_2 / R^(i+j),             (RT-3)
```

where `P_K` projects onto the scalar response basis. At `R=2`, the maximum of
RT-3 decreases from `0.263523` at K2 to `0.0862754` at K12. The raw complement
fraction can nevertheless be close to one: most of an individual interaction
may point into transverse response orbits even when its input-weighted size is
small.

Therefore the scalar tensor bound is only half of the required lemma. The
next estimate must resolve the complement by output shell, add the explicit
transverse `(2,1,1)` and oriented `(3,1,0)` orbit families, and prove a
cutoff-uniform weighted bound in a norm compatible with SLD-1. A plain L2
complement bound is not sufficient to control the derivative-weighted local
quartic expression.

## Graded transverse closure

The scalar response sequence does not contain every direction observed in the
trajectory extremizer. `LocalSldResponseBasis` therefore forms one common
graded orthonormal basis from the scalar responses, the transverse `(2,1,1)`
orbit, and both oriented `(3,1,0)` orbits. Lower-degree candidates are inserted
before higher response orders, preventing Gram--Schmidt from giving a
low-degree orbit a cutoff-dependent high-degree weight.

For several directions of the same analytic degree, an absolute coefficient
sum would depend on the arbitrary orthonormal basis. The augmented tensor
therefore also reports the degree-block invariant

```text
C_block(K) = max_(i,j) 1/R^(d_i+d_j)
  sum_d r^d [sum_(m:d_m=d) |T_mij|^2]^(1/2).                (RT-4)
```

A deterministic bilinear-closure step finds the largest weighted product
orthogonal to the current space, inserts its raw product at degree
`d_i+d_j+1`, and reorthonormalizes the complete basis in graded order. This is
a diagnostic construction of the missing response tree, not an assumption
that a finite number of directions closes the PDE.

With all three explicit orbits and sixteen closure extensions at `R=2`,
`r=23/20`, the current ledger is

```text
K   basis   C_block(K)   shell complement   combined finite bound
3    23     1.044672519      0.126203590          1.170876109
4    24     1.048938752      0.257011570          1.305950322
5    25     1.048938752      0.271164813          1.320103566
```

The combined column is the sum of separately maximized projected and
shell-complement terms, so it is a valid but non-sharp finite bound. The K4
and K5 projected constants agree exactly at serialized precision. The data
suggests the next stronger target:

```text
RT-G: the complete graded response/orbit tree has a cutoff-uniform
      two-radius bilinear norm from R=2 to r=23/20.          (RT-G)
```

RT-G still needs an analytic dimension/support count for every degree and a
summable estimate for the indefinitely generated transverse tree. Sixteen
greedy extensions and cutoffs through K5 do not replace either proof.

There is also a concrete obstruction to a support-only proof. Response degree
is a recursion label, not a strict Fourier convolution filtration. At K5 the
ledger finds output degree excess

```text
d_out - (d_left+d_right+1) = 2
```

with coefficient magnitude `0.409773`; one instance couples the axis and the
second transverse closure direction to `response-order-5`. Thus RT-G cannot
follow merely by declaring the tensor block-triangular. Its proof must
establish quantitative off-diagonal coefficient decay in addition to physical
Fourier-shell support. The certificate records `graded_support_closed=false`
so this rejected shortcut cannot be used silently.

## Reproduction

```bash
./build/navier_stokes_lab local-sld-response-tensor \
  --cutoff 12 --depth 13 \
  --input-radius 2 --output-radius 1.15 \
  --tolerance 1e-14 --threads 12 \
  --certificate proof/l4/analysis/shifted-local-density/response-tensor/R200-r115/K12.json
```

Primary artifacts:

- [`../../analysis/shifted-local-density/response-tensor/R200-r115/K12.json`](../../analysis/shifted-local-density/response-tensor/R200-r115/K12.json)
- [`../../analysis/shifted-local-density/response-tensor/R200-r125/K12.json`](../../analysis/shifted-local-density/response-tensor/R200-r125/K12.json)
- [`../../analysis/shifted-local-density/response-tensor/augmented-graded-R200-r115/K3-closure16.json`](../../analysis/shifted-local-density/response-tensor/augmented-graded-R200-r115/K3-closure16.json)
- [`../../analysis/shifted-local-density/response-tensor/augmented-graded-R200-r115/K4-closure16.json`](../../analysis/shifted-local-density/response-tensor/augmented-graded-R200-r115/K4-closure16.json)
- [`../../analysis/shifted-local-density/response-tensor/augmented-graded-R200-r115/K5-closure16.json`](../../analysis/shifted-local-density/response-tensor/augmented-graded-R200-r115/K5-closure16.json)

The direct-triad tensor, Gram errors below `5.7e-18`, exact norm
reconstruction, and the analytic axis constant are checked independently.
