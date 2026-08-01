# Equal-length orthogonal local triads

This directory isolates and controls the dominant squared-length signature
observed in the broad-spread local adversary. The result is an exact partial
lemma for one infinite family of local triads; it is not a proof of the whole
L4 lemma.

## Geometry

Consider ordered integer waves satisfying

```text
|p| = |q|,       p dot q = 0,       k = p+q.
```

Their sorted squared-length signature is `(r,r,2r)`. In a Fourier cube with
component cutoff `K`, fix a nonzero `p` and eliminate `q_z` using `p dot q=0`
after choosing a nonzero component of `p`. Substitution into `|q|^2=|p|^2`
gives a positive-definite quadratic equation in the two remaining coordinates.
For every fixed value of one coordinate there are at most two real values of
the other. Therefore

```text
degree(p) <= 2(2K+1).                                      (OT-1)
```

For a fixed target `k`, the identities imply

```text
2 p dot k = |k|^2,        |p|^2 = |k|^2/2.
```

The same plane/sphere argument proves the identical target-degree bound. The
bound is elementary and cutoff-uniform in form; the program enumeration is a
regression check, not its proof.

## Trilinear multiplicity bound

Let `E_R` be the Fourier energy in the two adjacent shells containing a
signature `(R^2,R^2,2R^2)`. For its interaction hypergraph, OT-1 gives maximum
degree `D_R <= C R`. For arbitrary nonnegative coefficient sequences,
Cauchy--Schwarz gives

```text
sum_edges a_p b_q c_k
 <= (sum_edges a_p^2 b_q^2)^(1/2)
    (sum_edges c_k^2)^(1/2)
 <= sqrt(D_R) ||a||_2 ||b||_2 ||c||_2.                    (OT-2)
```

The complete-triad frequency-spread lemma LS-1 contributes two derivatives,
and the advective coefficient contributes one. Combining OT-1 and OT-2 yields

```text
|V_orthogonal,R| <= C R^(7/2) E_R^(3/2).                  (OT-3)
```

The palinstrophy on these modes satisfies `P_R >= R^4 E_R`. Since
`E_R <= E(0)`,

```text
|V_orthogonal,R|
 <= C R^(-1/2) E(0)^(1/2) P_R.
```

For all sufficiently large `R`, with threshold depending only on `E(0)`,
viscosity, and the explicit finite constants, this is absorbed by `nu P_R/2`.
Only finitely many lower shells remain, and their contribution is bounded by
a finite-dimensional constant depending on the same data. Thus the entire
equal-length orthogonal signature family is not an independent blow-up
obstruction.

## Certificate

`OrthogonalTriadGeometry` stores the exact exponent calculation and enumerates
the lattice degrees as a regression. Through cutoff eight it reports

```text
K   ordered pairs   max input degree   max target degree   2(2K+1)
1       24                 4                   2                6
4      432                 4                   6               18
8     1728                12                   6               34
```

The maximum observed input/bound and target/bound ratios are `2/3` and `3/7`.
The exact scaling certificate gives transfer exponent `7/2` and
transfer-to-viscosity exponent `-1/2`. The machine-readable artifact is
[`orthogonal-triad-certificate.json`](orthogonal-triad-certificate.json).

Reproduce it with

```bash
./build/navier_stokes_lab orthogonal-triad-certificate \
  --max-cutoff 8 \
  --certificate proof/l4/lemmas/local-orthogonal-triads/orthogonal-triad-certificate.json
```

## Remaining gap

The K3 broad extremizer places `91.27%` of its signed local transfer in this
controlled family. That percentage is experimental and cannot discard the
remaining signatures. Nearby broad-spread length triples do not satisfy the
two exact quadratic constraints used in OT-1. The next task is to determine
whether they can be covered by a bounded number of thin signature classes with
the same subcritical degree estimate, or whether a genuinely dense local
family restores the critical `R^(9/2)` bound.

More precisely, if the coherent local interaction degree scales as `R^d`, the
same Cauchy argument gives frequency power `3+d/2`. The viscosity threshold is
`d=2`: the unrestricted local graph has the supercritical value `d=3`, while
the orthogonal family has the subcritical value `d=1`. The next target can
therefore be stated without ambiguity: prove an effective coherent degree
strictly below `R^2` for the unresolved broad-spread transfer, or find a
different cancellation before summing its dense signatures.
