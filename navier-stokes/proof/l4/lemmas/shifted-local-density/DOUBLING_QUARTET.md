# Doubling-family quartet closure

This note proves a cutoff-independent estimate for the complete dominant
`(m,m,2m)` local quartet. It closes the doubling-family block at the LQC-3
target scale. It does not close the signature remainder, the mixed block, or
the full local SLD lemma.

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

## Global sum of structural entries

A block whose low squared radius is `m` is supported only at squared radii
`m` and `2m`. Two outer, commutator, or nested structural blocks can therefore
pair only if their low radii differ by at most a factor two. In the
factor-four squared-radius ledger this is at most one neighboring shell.

Let `E_near,j` contain that fixed-width neighborhood. The one-shell result
and bounded overlap give

```text
sum_j R_j^5 E_near,j^2
    <= C Z^(3/2) P^(1/2)
    <= C Z^(5/4) P^(3/4).                            (DQ-3)
```

For the first inequality, write `z_j=R_j^2 E_near,j` and
`p_j=R_j^4 E_near,j`. Then
`R_j^5 E_near,j^2=z_j^(3/2)p_j^(1/2)` and use
`z_j<=sum z_i`, `p_j<=sum p_i`, plus bounded overlap. The second inequality
is `(Z/P)^(1/4)<=1`, which follows from the mean-zero torus spectral gap.

## Rejected intermediate absorption

DQ-2 does not by itself sum the two global normalization entries. The tempting
sequence inequality is

```text
S = sum_j R_j^(7/2) E_j^(3/2),
T = sum_j R_j^(11/2) E_j^(3/2),
P = sum_j R_j^4 E_j,
M = sum_j R_j^5 E_j^2,

S T <= P M.                                             (DQ-4)
```

DQ-4 is false. Take one unit-energy shell at radius one and a second shell at
radius `L` with energy `L^(-11/4)`. The high-shell powers are

```text
S_high = L^(-5/8),    T_high = L^(11/8),
P_high = L^(5/4),     M_high = L^(-1/2).
```

The low shell controls `S` and `M`, while the high shell controls `T` and `P`.
Consequently

```text
S T/(P M) grows like L^(1/8).                          (DQ-5)
```

The C++ two-shell ledger independently reaches ratio `2.99973` within its K12
scan and emits the maximizing tested configuration. DQ-5 is the analytic
counterexample; the finite scan only verifies the implementation.

Thus the normalization term cannot first be absorbed into the structural
shell sum `M`. This rejects that intermediate proof route, but it does not
reject the target estimate itself.

## Direct normalization estimate

The normalization sums can instead be sent directly to the LQC-3 target.
The incidence bound gives

```text
|S_j| <= C R_j^(7/2) E_near,j^(3/2)
|T_j| <= C R_j^(11/2) E_near,j^(3/2),
```

where `T=<A B_d,A u>`. In terms of `z_j` and `p_j`,

```text
R_j^(7/2) E_j^(3/2)  = z_j^(5/4) p_j^(1/4),
R_j^(11/2) E_j^(3/2) = z_j^(1/4) p_j^(5/4).
```

Since each `z_j<=Z` and `p_j<=P`, bounded overlap yields

```text
|S| <= C Z^(5/4) P^(1/4),
|T| <= C Z^(1/4) P^(5/4).                            (DQ-6)
```

Therefore both global normalization entries close:

```text
S^2/Z <= C Z^(3/2)P^(1/2),
|S T|/P <= C Z^(3/2)P^(1/2)
         <= C Z^(5/4)P^(3/4).                        (DQ-7)
```

Together, DQ-3 and DQ-7 prove

```text
|K_d+G_d| <= C Z^(5/4)P^(3/4),                       (DQ-8)
```

with `C` independent of the Galerkin cutoff. This is stronger than the
doubling-block instance of LQC-3 up to its fixed positive initial-data factor.
The exact-rational engine checks both exponent routes independently. The
failed DQ-4 absorption screen remains in the certificate so that this invalid
shortcut is not accidentally reintroduced.

## Exact projected-square diagnostic

`LocalSldProjectedSquare` also verifies, for any fixed triad selection,

```text
-<B,A B> + S^2/(2Z) + 3S T/(2P)
 = -||A^(1/2)(B-c A u)||_2^2
   + S^2/(2Z) + c^2 H3,
c=3S/(4P).
```

This identity is useful for retaining sign, but DQ-5 is what removes the
otherwise non-closing `H3` remainder. The shell ledger reconstructs the
complete ordered matrix and the completed square with relative errors around
`10^-19` on the stored examples.

A 275-state two-scale screen over `L=2,...,12` and 25 response angles tests
the rejected energy profile `E_high=L^(-11/4)` directly. The largest positive
target ratio is `0.016067970718` at `L=2`; the per-scale positive maximum
decreases to `0.003136407508` at `L=12`. This finite scan is a regression and
falsification artifact, not the proof of DQ-8.

## Reproduction

```bash
./build/navier_stokes_lab doubling-quartet-certificate \
  --max-cutoff 12 \
  --certificate proof/l4/analysis/shifted-local-density/doubling-quartet/closed-family-K12.json

./build/navier_stokes_lab local-sld-doubling-scale-scan \
  --min-scale 2 --max-scale 12 \
  --angle-min -1.2 --angle-max 1.2 --angle-count 25 \
  --energy-decay-power 2.75 --threads 12 \
  --certificate proof/l4/analysis/shifted-local-density/doubling-quartet/two-scale-L2-L12-angle-scan.json
```

Artifact:

- [`../../analysis/shifted-local-density/doubling-quartet/closed-family-K12.json`](../../analysis/shifted-local-density/doubling-quartet/closed-family-K12.json)
- [`../../analysis/shifted-local-density/doubling-quartet/two-scale-L2-L12-angle-scan.json`](../../analysis/shifted-local-density/doubling-quartet/two-scale-L2-L12-angle-scan.json)
- [`../../analysis/shifted-local-density/doubling-quartet/two-scale-L12-Ehigh-Lm11over4-shell-matrix.json`](../../analysis/shifted-local-density/doubling-quartet/two-scale-L12-Ehigh-Lm11over4-shell-matrix.json)

The certificate reports `cutoff_independent_closed_family_bound=true` and,
deliberately, `full_local_lemma_proved=false`. The next proof block is the
closed signature remainder, followed by the mixed term.
