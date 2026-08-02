# Equal-low triple-family quartet closure

This note closes the complete squared-length family `(m,m,3m)` at the LQC-3
scale. It is a cutoff-independent sub-lemma. It does not close the sum over
all other projective shapes, the mixed block, or the full shifted local-density
lemma.

## Incidence geometry

Let `p+q=r` with

```text
|p|^2=|q|^2=m,   |r|^2=3m.
```

Then

```text
p dot q=m/2,   p dot r=3m/2.
```

Thus, after one input or the target is fixed, the other input lies on the
intersection of one integer plane and one sphere. Choose a nonzero coordinate
of the fixed vector and eliminate it with the plane equation. Fixing one of
the two remaining coordinates leaves a quadratic equation, hence at most two
solutions. The input and target degrees are therefore `O(R)` inside a cube of
radius `R`. This argument holds at every cutoff; the finite K12 enumeration is
only a regression test for the indexing implementation.

Target-wise Cauchy--Schwarz gives

```text
||B_3(v,w)||_2 <= C R^(3/2)||v||_2||w||_2.          (TQ-1)
```

## Quartet power

The two derivatives outside the two bilinear factors contribute `R^2`, so

```text
|K_3+G_3|_j <= C R_j^5 E_near,j^2.                 (TQ-2)
```

The LQC-3 shell target is `R_j^(11/2)E_near,j^2`. TQ-2 gains one half
derivative. The fixed ratio `sqrt(3)` makes the structural terms local in a
fixed number of neighboring dyadic shells. The same direct sequence estimates
used for the doubling family give

```text
S_3 <= C Z^(5/4)P^(1/4),
T_3 <= C Z^(1/4)P^(5/4).                            (TQ-3)
```

Consequently the complete `(m,m,3m)` family, including its normalization
entries, is bounded cutoff-independently at `Z^(5/4)P^(3/4)`.

## Reproduction

```bash
./build/navier_stokes_lab equal-low-quartet-certificate \
  --multiplier 3 --max-cutoff 12 \
  --certificate proof/l4/analysis/shifted-local-density/equal-low-quartet/triple-family-K12.json
```

Artifact:

- [`../../analysis/shifted-local-density/equal-low-quartet/triple-family-K12.json`](../../analysis/shifted-local-density/equal-low-quartet/triple-family-K12.json)

The certificate records `R^5` against the `R^(11/2)` target and keeps
`full_local_lemma_proved=false`.
