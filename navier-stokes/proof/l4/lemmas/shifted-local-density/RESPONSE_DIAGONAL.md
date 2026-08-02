# Cutoff-diagonal response majorant

This note removes a Galerkin-boundary artifact from the finite quadratic
response hierarchy and states the next explicit analytic target. It does not
prove the Clay problem or the local SLD lemma.

## Boundary-free diagonal

Let `b_0` be the normalized cyclic axis state and construct orthonormal
quadratic responses from

```text
raw b_n = sum_{i+j=n-1} B(b_i,b_j).
```

At Galerkin cutoff `K`, response order `n>K` is contaminated because its raw
convolution can require shell `K+1`. Fixed depth-16 comparisons at K3 and K4
therefore mix genuine cutoff convergence with a boundary-reflected tail.

The cutoff-diagonal ledger retains only orders `0,...,K`. For a normalized
reference state `u`, define

```text
a_n(K) = |<u,b_n>| / sqrt(E),
A_r(K) = sum_{n=0}^K r^n a_n(K),       r >= 1.              (RD-1)
```

The first omitted order `K+1` is still constructed and reported, but it is a
boundary diagnostic rather than part of `A_r(K)`.

## Exact sequence algebra

For nonnegative coefficients, the quadratic Cauchy product obeys

```text
sum_{n=1}^K r^n sum_{i+j=n-1} a_i a_j
    <= r A_r(K)^2.                                           (RD-2)
```

`LocalSldResponseDiagonal` evaluates both sides. RD-2 is exact elementary
algebra; it is not yet a Navier--Stokes estimate because the norms and
projections of `B(b_i,b_j)` still require a cutoff-uniform bound.

## Current numerical diagonal

K2 and K3 are optimized frozen-trajectory states. K4 is the six-step,
12-restart continuation from K3. K5 through K8 are exact zero-padding of that
same K4 state; they test the response construction, not fresh extremizers.

For `r=1.25`:

```text
K   safe orders   A_r(K)       weighted tail n>=3   |a_K|
2       3         1.120120280       0                1.01603e-2
3       4         1.142036207       1.88634e-2       9.65808e-3
4       5         1.148628944       2.51370e-2       2.58799e-3
5       6         1.151472355       2.79804e-2       9.31729e-4
6       7         1.151915110       2.84232e-2       1.16066e-4
7       8         1.158370374       3.48784e-2       1.35377e-3
8       9         1.159159126       3.56672e-2       1.32331e-4
```

The common safe coefficient difference is `1.67e-4` from K3 to the optimized
K4 state. It is at most `1.03e-21` for the common orders of the identical
zero-padded K4 state at K5 through K8. The weighted sums at K8 are
`1.10975119`, `1.15915913`, `1.23403534`, and `1.58091278` for radii
`1`, `1.25`, `1.5`, and `2`, respectively.

Artifacts:

- [`../../analysis/shifted-local-density/cyclic-response-diagonal-K2-K8-r100.json`](../../analysis/shifted-local-density/cyclic-response-diagonal-K2-K8-r100.json)
- [`../../analysis/shifted-local-density/cyclic-response-diagonal-K2-K8-r125.json`](../../analysis/shifted-local-density/cyclic-response-diagonal-K2-K8-r125.json)
- [`../../analysis/shifted-local-density/cyclic-response-diagonal-K2-K8-r150.json`](../../analysis/shifted-local-density/cyclic-response-diagonal-K2-K8-r150.json)
- [`../../analysis/shifted-local-density/cyclic-response-diagonal-K2-K8-r200.json`](../../analysis/shifted-local-density/cyclic-response-diagonal-K2-K8-r200.json)

## Required lemma

The finite data motivates, but does not establish, the following route:

```text
RD-L: for some r>1, the relevant response coefficients obey
      sup_K A_r(K) < infinity,

and the response basis satisfies a cutoff-independent bilinear majorant

      ||projection_m B(b_i,b_j)|| <= beta_{i,j,m}

whose r-weighted convolution norm is finite.
```

To imply the local SLD estimate, RD-L must be supplemented by a bound for the
orthogonal complement of the response/orbit space and must hold along every
Galerkin trajectory from the fixed smooth datum, not merely on optimized
finite states. Proving either requirement by assuming a uniform future
Sobolev norm would be circular.

The next machine task is a sparse interaction-tensor ledger for
`<b_m,B(b_i,b_j)>` on the boundary-free diagonal. The next analytic task is to
turn its exact symmetry and shell support into a cutoff-independent weighted
operator bound.
