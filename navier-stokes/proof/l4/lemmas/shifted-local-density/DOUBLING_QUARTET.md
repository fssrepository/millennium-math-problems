# Doubling-family quartet closure

This note gives a cutoff-independent one-shell estimate for the dominant
`(m,m,2m)` local quartet and records the exact cross-shell obstruction that
remains. It closes neither the complete doubling-family bracket nor the local
SLD lemma.

## Orthogonal incidence bound

The signature `(m,m,2m)` is equivalent to

```text
|p|=|q|,       p dot q=0,       k=p+q.
```

At frequency `R`, fixing an input or target leaves integer points on the
intersection of one plane and one sphere. Choosing one free coordinate and
solving the remaining quadratic gives at most `2(2R+1)` possibilities. The
Leray projection has norm one, so target-wise Cauchy--Schwarz yields

```text
||B_d(v,w)||_2 <= C R^(3/2) ||v||_2 ||w||_2.         (DQ-1)
```

The `R^(3/2)` consists of one derivative and the square root of the `O(R)`
incidence degree. The finite lattice implementation checks the elementary
degree bound through K12; the preceding plane--sphere argument is the
cutoff-independent proof, not the enumeration.

## Complete one-shell power

Let `E` be the energy in one dyadic neighborhood. DQ-1 gives

```text
||B_d(u,u)||_2 <= C R^(3/2) E.
```

Every entry in the closed `K_d+G_d` bracket then has the same scale:

```text
entry                                      frequency scale
<B_d,A B_d>                                R^5 E^2
<B_d,B_d(u,A u)>                           R^5 E^2
<A u,B_d(-B_d,u)>                          R^5 E^2
S_d^2/Z                                    R^5 E^2
S_d <A B_d,A u>/P                          R^5 E^2
```

The LQC-3 target has one-shell scale `R^(11/2)E^2`. Thus the orthogonal
doubling family has a genuine half-derivative gain:

```text
R^5 E^2 / [R^(11/2)E^2] = R^(-1/2).        (DQ-2)
```

All exponents are evaluated as exact rationals by
`DoublingQuartetClosure`; the self-test requires `3/2`, `5`, `11/2`, and
`-1/2` exactly.

## Rejected absolute cross-shell sum

DQ-2 does not by itself sum the two global normalization entries. The tempting
sequence inequality is

```text
S = sum_j R_j^(7/2) E_j^(3/2),
T = sum_j R_j^(11/2) E_j^(3/2),
P = sum_j R_j^4 E_j,
M = sum_j R_j^5 E_j^2,

S T <= P M.                                             (DQ-3)
```

DQ-3 is false. Take one unit-energy shell at radius one and a second shell at
radius `L` with energy `L^(-11/4)`. The high-shell powers are

```text
S_high = L^(-5/8),    T_high = L^(11/8),
P_high = L^(5/4),     M_high = L^(-1/2).
```

The low shell controls `S` and `M`, while the high shell controls `T` and `P`.
Consequently

```text
S T/(P M) grows like L^(1/8).                          (DQ-4)
```

The C++ two-shell ledger independently reaches ratio `2.99973` within its K12
scan and emits the maximizing tested configuration. DQ-4 is the analytic
counterexample; the finite scan only verifies the implementation.

Therefore the dominant family cannot be closed by taking absolute values of
the palinstrophy-normalization term separately. The next lemma must combine
its signed cross-shell contribution with the commutator, outer-square, and
advecting entries before summation. This is exactly the cancellation retained
by `K_d+G_d` and by the full block factorization.

## Reproduction

```bash
./build/navier_stokes_lab doubling-quartet-certificate \
  --max-cutoff 12 \
  --certificate proof/l4/analysis/shifted-local-density/doubling-quartet/closed-family-K12.json
```

Artifact:

- [`../../analysis/shifted-local-density/doubling-quartet/closed-family-K12.json`](../../analysis/shifted-local-density/doubling-quartet/closed-family-K12.json)

The certificate deliberately reports `full_local_lemma_proved=false` and
`cutoff_independent_closed_family_bound=false`.
