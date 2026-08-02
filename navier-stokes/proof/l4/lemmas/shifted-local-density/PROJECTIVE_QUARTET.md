# Fixed projective-ray quartet closure

This note extends the fixed-angle argument to every primitive squared-length
shape `(a,b,c)`. For each fixed shape it closes all integer dilates
`(am,bm,cm)`, `m>=1`, with a cutoff-independent constant that may depend on
the shape. It does not prove a uniform or summable estimate over the growing
set of primitive shapes.

## Plane--sphere incidence

Sort `a<=b<=c`, assume `gcd(a,b,c)=1`, and consider any ordered role
assignment

```text
|p|^2=am,   |q|^2=bm,   |p+q|^2=cm.
```

The real triangle condition is

```text
(c-a-b)^2 <= 4ab.                                   (PQ-1)
```

For fixed `p`, the scale is fixed by `m=|p|^2/a` and

```text
p dot q=(c-a-b)m/2.                                 (PQ-2)
```

Therefore `q` lies on one plane and one sphere. Eliminating a coordinate for
which `p` is nonzero and then fixing one free coordinate leaves at most two
roots. The input degree is at most `C K`. For fixed target `r=p+q`, the scale
is `m=|r|^2/c` and

```text
r dot p=(c+a-b)m/2,                                 (PQ-3)
```

so the same argument bounds the target degree. Summing the at most six role
assignments changes only the constant. This proves an `O(K)` degree bound for
the entire scale ray, not merely for one exact signature.

## Fixed-ray quartet estimate

PQ-2--PQ-3 give

```text
||B_(a,b,c)(v,w)||_2
  <= C_(a,b,c) K^(3/2)||v||_2||w||_2.               (PQ-4)
```

The closed quartet therefore has power

```text
2 + 2(3/2) = 5,
```

while the LQC-3 target has power `11/2`. For a fixed shape, the radius ratios
are fixed, so only finitely many dyadic-shell offsets occur. The direct `S`
and `T` sequence bounds then close the normalization terms as well. Hence

```text
|K_(a,b,c)+G_(a,b,c)|
  <= C_(a,b,c) Z^(5/4)P^(3/4)                       (PQ-5)
```

with a constant independent of the Galerkin cutoff.

## Square function and synthesis obstruction

The remainder contains all primitive shapes simultaneously, and the quadratic
quartet expansion also contains cross-ray terms. Applying PQ-5 shape by shape
and summing absolute values can reproduce the dense `R^(3/2)` loss.

There is one further cutoff-independent gain. Ordered input pairs partition
disjointly by primitive ray, and the plane--sphere degree constant is uniform
over local shapes. Target-wise Cauchy--Schwarz followed by summation over the
disjoint pair sets gives the exact square-function estimate

```text
sum_sigma ||B_sigma(u,u)||_2^2 <= C K^3 E^2.        (PQ-6)
```

Thus the projective square function itself retains the required `K^(3/2)`
bilinear scale without a shape-count loss. Its formal diagonal quartet power
would be `K^5`. A standalone coherent synthesis estimate, however, is not
available: it is false with a cutoff-independent constant.

For every sufficiently large integer `K`, set

```text
r=(0,0,K),  q=(0,a,b),  p=(0,-a,K-b),
a/K in [1/2,3/5],       b/K in [1/4,2/5].           (PQ-7)
```

Choose `u_q=e_x`, `u_p=i(0,-p_z,p_y)/|p|`, and impose conjugate symmetry at
the negative modes. All three radii are local, `r` has the unique largest
squared radius, and the primitive projective shape determines `(a,b)`.
There are at least `K^2/160` distinct primitive rays. At the common target
`r`, all their advection contributions point in the same positive `e_x`
direction and have coefficients in `[K/2,K]`. Consequently

```text
||sum_sigma B_sigma(u,u)_r||_2^2
-------------------------------- >= K^2/640.        (PQ-8)
 sum_sigma ||B_sigma(u,u)_r||_2^2
```

This proves that no cutoff-independent bound can synthesize the projective
square function by itself. The same construction also satisfies
`<Au,B(u,u)>=0` exactly: `Au` is supported in the `yz` polarization, whereas
the coherent target is in `e_x`. Hence its full stretching and its power-one
product both vanish. The obstruction rejects only the stronger standalone
synthesis route; it does not reject the remainder lemma.

The unresolved statement must therefore couple coherent synthesis to the
stretching direction tested by `Au`, or exploit an equivalent signed
cross-ray cancellation compatible with

```text
|(K_rem+G_rem)S_full| <= C Z^2P^2.                  (PQ-9)
```

The projective ledger groups exact signatures by their gcd. On the K7
power-one winner it reduces 8105 exact signatures to 7322 primitive shapes,
but one ray contributes `95.42%` of the absolute total. After removing the
proved doubling and triple families and the exact `(1,2,3)` signature, the K5
tail has only `3.027` effective projective shapes. These are finite diagnostics,
not the missing sum theorem.

## Reproduction

```bash
./build/navier_stokes_lab projective-quartet-certificate \
  --signature 2,3,5 --max-cutoff 12 \
  --certificate proof/l4/analysis/shifted-local-density/projective-quartet/ray-2-3-5-K12.json

./build/navier_stokes_lab projective-square-function-certificate \
  --certificate proof/l4/analysis/shifted-local-density/projective-quartet/square-function.json

./build/navier_stokes_lab projective-fan-certificate \
  --max-cutoff 1024 \
  --certificate proof/l4/analysis/shifted-local-density/projective-quartet/coherent-fan-obstruction-K1024.json

./build/navier_stokes_lab local-sld-projective-fan-scan \
  --min-cutoff 3 --max-cutoff 24 --threads 12 \
  --certificate proof/l4/analysis/shifted-local-density/projective-quartet/coherent-fan-unique-rays-K3-K24.json \
  --state-dir proof/l4/states/projective-coherent-fan/unique-rays-K3-K24
```

Artifacts:

- [`../../analysis/shifted-local-density/projective-quartet/ray-2-3-5-K12.json`](../../analysis/shifted-local-density/projective-quartet/ray-2-3-5-K12.json)
- [`../../analysis/shifted-local-density/projective-quartet/ray-3-3-8-K12.json`](../../analysis/shifted-local-density/projective-quartet/ray-3-3-8-K12.json)
- [`../../analysis/shifted-local-density/projective-quartet/ray-1-3-4-K12.json`](../../analysis/shifted-local-density/projective-quartet/ray-1-3-4-K12.json)
- [`../../analysis/shifted-local-density/projective-quartet/square-function.json`](../../analysis/shifted-local-density/projective-quartet/square-function.json)
- [`../../analysis/shifted-local-density/projective-quartet/coherent-fan-obstruction-K1024.json`](../../analysis/shifted-local-density/projective-quartet/coherent-fan-obstruction-K1024.json)
- [`../../analysis/shifted-local-density/projective-quartet/coherent-fan-unique-rays-K3-K24.json`](../../analysis/shifted-local-density/projective-quartet/coherent-fan-unique-rays-K3-K24.json)
- [`../../adversary/shifted-local-density/projective-coherence/K1-K5.json`](../../adversary/shifted-local-density/projective-coherence/K1-K5.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K7-power-one-projective.json`](../../analysis/shifted-local-density/remainder-quartet/K7-power-one-projective.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K5-power-one-tail-projective.json`](../../analysis/shifted-local-density/remainder-quartet/K5-power-one-tail-projective.json)

The analytic plane--sphere argument is the proof of each fixed-ray statement;
the K12 counts verify only the corresponding discrete implementation.
