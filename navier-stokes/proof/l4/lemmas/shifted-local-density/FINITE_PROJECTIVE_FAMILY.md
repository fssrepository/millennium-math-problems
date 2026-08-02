# Fixed finite projective-family closure

This note closes every quartet term whose projective rays all belong to one
fixed finite family. It includes unequal-ray cross terms inside that family.
It does not control coupling to a growing family of rays and does not prove
the full shifted local-density lemma.

## Statement

Let `F` be any nonempty finite set of primitive feasible squared-length
signatures and write

```text
B_F(v,w) = sum_(sigma in F) B_sigma(v,w),
b_F      = B_F(u,u),
c_F      = B_F(u,Au),
S_F      = <Au,b_F>,
T_F      = <A b_F,Au>.
```

The quartet internal to `F` is

```text
J_F = -<b_F,A b_F>+<b_F,c_F>-<Au,B_F(b_F,u)>
      +S_F^2/(2Z)+3S_F T_F/(2P).                  (FPF-1)
```

Expanding `B_F=sum B_sigma` shows that FPF-1 is exactly the sum of every
self term PCQ-1 and every unequal-ray cross term PCQ-3 with both rays in
`F`. No absolute sum over ray pairs is introduced.

For every fixed `F`, there is a constant `C_F`, independent of the Galerkin
cutoff, such that

```text
|J_F| <= C_F Z^(5/4) P^(3/4).                     (FPF-2)
```

## Proof

For one primitive ray, fixing an input or target reduces every role assignment
to a plane--sphere intersection. Its degree is `O_sigma(K)`. A finite union
therefore has

```text
degree(B_F) <= sum_(sigma in F) C_sigma K = C_F K. (FPF-3)
```

The standard bounded-degree convolution estimate gives

```text
||B_F(v,w)||_2 <= C_F K^(3/2)||v||_2||w||_2.      (FPF-4)
```

Use FPF-4 directly on the unexpanded operator FPF-1. The structural terms
have shell-frequency power

```text
2+2(3/2)=5,
```

whereas the required LQC-3 power is `11/2`. Thus the structural block has a
half-derivative gain. Every fixed ray has only finitely many radial shell
offsets. Their finite union has the same property, and finite
Cauchy--Schwarz applied to the direct sequence estimates gives

```text
|S_F| <= C_F Z^(5/4)P^(1/4),
|T_F| <= C_F Z^(1/4)P^(5/4).                      (FPF-5)
```

Consequently `S_F^2/Z` and `S_F T_F/P` also have power five and satisfy
FPF-2. This proves the complete internal self+cross bound. The constant may
grow with `F`; finiteness is essential.

## Remaining boundary

For a fixed core `F` and its growing complement `T`, the exact projective
quartet still contains

```text
J_all = J_F + J_FT + J_T.                          (FPF-6)
```

FPF-2 closes only `J_F`. The open problem is a signed, stretching-aware
estimate for the core--tail term `J_FT` and the tail--tail term `J_T` after
multiplication by `S_full/(Z^2P^2)`. Allowing `F` to grow with the cutoff
would merely rename this missing estimate.

## Canonical height cores and exact boundary ledger

For an integer `H`, define the canonical core

```text
F_H = {(a,b,c): 1<=a<=b<=c<=H, gcd(a,b,c)=1,
                   (c-a-b)^2<=4ab}.                (FPF-7)
```

`F_H` is fixed and finite, so FPF-2 applies. `ProjectiveCoreFamily`
enumerates exactly this set. The exact boundary ledger evaluates all five
components of FPF-1 separately for `J_F`, `J_FT`, and `J_T`, using 12-thread
projective interaction sums. On the original K7 power-one winner, the
13-ray empirical core gives

```text
selected       fixed core       core--tail       tail--tail
0.0063647724   0.0063635050     3.1639e-7        9.5097e-7.
```

The bracket reconstruction error is `6.42e-19`. This localizes the open
quantity exactly; it is not a tail estimate.

The height-core scan initially showed rapid decay on old winners, but the
exact-gradient objective

```text
|(J_FT+J_T) S_full|^2/(Z^4P^4)                    (FPF-8)
```

finds states adapted to the boundary. Its gradient agrees with central
differences to `8.94e-12`. On separately optimized K8 states, the square
roots of FPF-8 are

```text
H                     8          16          32          64
open power one   0.0010076   0.0008591   0.0006330   0.0003079
```

The fitted finite lower-branch height exponent is `-0.557`. Tail--tail is
larger at H=8 and H=16; by H=64 the core--tail and tail--tail pieces are of
comparable size. The observed decay is encouraging but is neither an upper
bound nor an asymptotic theorem: each row used only one finite-cutoff local
optimizer, and the optimizer visibly moves toward higher primitive shapes as
`H` increases.

## Machine certificate

The C++ certificate canonicalizes and rejects duplicate or nonprimitive
signatures, applies the exact plane--sphere argument member by member, sums
their theoretical incidence bounds, and audits the discrete degrees through
a requested cutoff. The audit is a consistency check; FPF-3--FPF-5 are the
cutoff-independent proof.

```bash
./build/navier_stokes_lab finite-projective-family-certificate \
  --max-cutoff 8 \
  --signature 1,1,1 --signature 1,1,3 --signature 1,2,3 \
  --signature 1,3,3 --signature 1,3,4 --signature 2,3,5 \
  --signature 2,5,5 --signature 3,3,4 --signature 3,3,8 \
  --signature 3,4,11 --signature 3,5,6 --signature 3,8,11 \
  --signature 5,5,6 \
  --certificate proof/l4/certificates/shifted-local-density/projective/finite-dominant-core-K8.json

./build/navier_stokes_lab finite-projective-family-certificate \
  --max-cutoff 8 --maximum-height 8 \
  --certificate proof/l4/certificates/shifted-local-density/projective/finite-height-8-core-K8.json

./build/navier_stokes_lab local-sld-projective-open-power-replay \
  --height 8 --state proof/l4/states/local-projective-open-power/height-8-K8-refined/K8.tsv \
  --height 16 --state proof/l4/states/local-projective-open-power/height-16-K8-refined/K8.tsv \
  --height 32 --state proof/l4/states/local-projective-open-power/height-32-K8-warm/K8.tsv \
  --height 64 --state proof/l4/states/local-projective-open-power/height-64-K8-warm/K8.tsv \
  --threads 12 --exclude-triple-family \
  --certificate proof/l4/analysis/shifted-local-density/remainder-quartet/K8-optimized-open-power-height-replay.json
```

The listed 13-ray core contains the leading rays found by the power-one and
cross-attribution ledgers. Its certificate records internal self+cross
closure with frequency gain `-1/2`, while explicitly recording the
core--tail, tail--tail, uniform-growing-family, full-lemma, and Clay claims as
unproved.

Artifacts:

- [`../../certificates/shifted-local-density/projective/finite-dominant-core-K8.json`](../../certificates/shifted-local-density/projective/finite-dominant-core-K8.json)
- [`../../certificates/shifted-local-density/projective/finite-height-8-core-K8.json`](../../certificates/shifted-local-density/projective/finite-height-8-core-K8.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K7-power-one-dominant-core-tail.json`](../../analysis/shifted-local-density/remainder-quartet/K7-power-one-dominant-core-tail.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K8-optimized-open-power-height-replay.json`](../../analysis/shifted-local-density/remainder-quartet/K8-optimized-open-power-height-replay.json)
