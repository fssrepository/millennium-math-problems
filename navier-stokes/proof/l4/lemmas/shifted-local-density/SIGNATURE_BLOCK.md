# Exact dominant-signature decomposition

The direct SLD adversary identifies the squared-frequency signature
`(1,1,2)` as the dominant low-shell interaction. This file records an exact
operator split; it does not infer a proof from the observed dominance.

Let `L` be the local bilinear advection operator and split it as

```text
L = L_d + L_r,

L_d: sorted (|p|^2,|q|^2,|k|^2)=(1,1,2),
L_r: every other local signature.
```

The code also supports the scale family `(m,m,2m)`, for all positive squared
lengths `m`, and its exact local complement. Write

```text
B_d=L_d(u,u),  B_r=L_r(u,u),
S_d=<Au,B_d>,  S_r=<Au,B_r>,
T_d=<AB_d,Au>, T_r=<AB_r,Au>.
```

The closed dominant and remainder brackets are obtained from the same `K+G`
formula with `L,B,S,T` replaced consistently by their `d` or `r` versions.
The independently evaluated mixed bracket is

```text
F_x = -<B_d,AB_r>-<B_r,AB_d>
      +<B_d,L_r(u,Au)>+<B_r,L_d(u,Au)>
      +S_d S_r/Z
      +3(S_d T_r+S_r T_d)/(2P)
      -<Au,L_d(B_r,u)>-<Au,L_r(B_d,u)>.
```

Consequently

```text
K+G = F_d+F_r+F_x.                                  (LSB-1)
```

Multiplication by the common signed factor

```text
4S^3 ZP/[k0(S^4Z^2P+B0Z^3P^4)]
```

gives an exact three-way decomposition of `R_local` without changing the
normalization between blocks.

For the K8 direct-SLD winner, the exact `(1,1,2)` split is

```text
full R_local                 7.95364960092e-4
closed dominant block        7.70181961863e-4
closed remainder             1.21606302797e-5
mixed block                  1.30223679496e-5
dominant absolute fraction   0.968337807808
reconstruction error         2.20e-18
```

Selecting every `(m,m,2m)` signature changes the dominant value by only
`4e-15` on this state. This says that higher members of the scale family are
inactive on the observed branch; it does not permit their omission in a
general proof.

The next conventional estimate must treat two parts:

1. a scale-uniform signed bound for the complete `(m,m,2m)` family;
2. a summable bound for `F_r+F_x` which preserves the mixed signs in LSB-1.

[`DOUBLING_QUARTET.md`](DOUBLING_QUARTET.md) now proves the
cutoff-independent one-shell power for item 1. Orthogonal incidence gives
`|K_d+G_d|_j <= C R_j^5 E_{near,j}^2`, a half-derivative improvement over the
required `R_j^(11/2)E_{near,j}^2` scale. A two-scale construction also proves
that separately summing the absolute palinstrophy-normalization term is false:
the naive sequence ratio grows like `L^(1/8)`. Item 1 is therefore reduced to
a signed cross-shell cancellation inside the already combined `K_d+G_d`, not
an unresolved one-shell counting estimate.

Bounding `F_d`, `F_r`, and `F_x` by unrelated raw norms would discard the
decomposition's purpose and can restore the previously identified
`(P/Z)^(3/4)` loss.

## Common-normalization block objectives

The block engine now differentiates the three actual contributions, rather
than replacing the full stretching `S` by a block stretching. Define

```text
c_j = F_j E^(1/4)/(Z^(7/4)P),
x   = S/(E^(1/4)Z^(1/4)P),
phi(x)=4x^3/(1+x^4),
R_j=c_j phi(x),                    j=d,r,x.
```

Here `S` and `phi` are always computed from the full local operator. Reverse
mode uses

```text
grad R_j = phi grad c_j + c_j phi'(x) grad x,
phi'(x)  = 4x^2(3-x^4)/(1+x^4)^2,
c_x      = c_full-c_d-c_r.
```

Central differences validate the dominant, remainder, mixed, and reconstructed
full gradients at relative errors between `7.15e-12` and `1.29e-11`. The
three ratios reconstruct the full ratio to `1.97e-19` in the regression.

Independent 12-start projected L-BFGS scans give:

```text
block                    K1             K2             last measured
doubling family       7.63980e-4     7.72930e-4     7.72930e-4 (K6)
local remainder       1.54701e-4     2.20222e-4     2.20683e-4 (K6)
mixed                 1.21920e-4     1.34752e-4     1.34936e-4 (K6)
```

These maxima occur at different states and must not be added as if they were
one observed solution. They show that the remainder and mixed terms are not
uniformly negligible merely because they are small on the K8 full-objective
winner. All three observed branches are low-mode and nearly cutoff stable;
that is evidence for the shape of a bound, not a proof of one.

## Frozen trajectory peak

The decomposition can evolve an optimized initial state and retain the
initial `k0,B0` in the common denominator. Warm continuation from K2 moves the
K3 `[0,0.5]` maximum to the stronger branch at `t=0.298`, where it gives

```text
full frozen R_local             8.53498799310e-4
closed doubling family          8.24077287430e-4
closed local remainder          1.51420260355e-5
mixed block                     1.42794858448e-5
dominant absolute fraction      0.965528350006
reconstruction error            6.20e-20
```

Thus the same family remains dominant after evolution with the correct frozen
normalization. The two corrections remain comparable and cannot be dropped.
The winning initial state has `99.969%` of its energy in the first hard shell,
while the `(1,1,2)` signature carries `99.627%` of coherent local transfer.
This routes the main estimate through the complete `(m,m,2m)` family and a
small, explicitly structured response hierarchy rather than an arbitrary
high-dimensional state.

The hierarchy has now been made explicit: 16 Gram--Schmidt response orders,
one transverse `(2,1,1)` polarization, and two oriented `(3,1,0)` cyclic
orbits capture `99.9998628%` of the current K3 winner's energy. Their projected
state retains `99.96181%` of the frozen trajectory objective at the same peak
time. This supports an orbit-wise summability lemma, but finite K3 projection
accuracy does not establish cutoff-uniform control of `F_d+F_r+F_x`.

## Rejected sign shortcut

The stronger shortcut `K+G<=0` is false. The signed-closure gradient has
central-difference error `1.57e-11`, and a 12-start search finds

```text
full local c                       0.0740187069851
doubling-family c                  0.0729937
doubling-complement c              0.0740187
```

The full counterexample has essentially six active Fourier modes and
`S=8.78e-9` after continuation. Thus it decisively rejects a sign theorem for
the bracket, while its `S^3`-weighted SLD contribution is negligible. The
viable target must retain the joint `c_j phi(x)` structure.

Reproduce the scale-family certificate with:

```bash
./build/navier_stokes_lab local-sld-block \
  --state proof/l4/states/local-sld-ratio/K8/K8.tsv \
  --doubling-family --threads 12 \
  --certificate proof/l4/analysis/shifted-local-density/direct-local-sld-K8-doubling-family.json
```

The block adversaries are reproducible with `local-closure-adversary` using
`--objective block-ratio` or `--objective mixed-ratio` and
`--selection doubling-family` (or `doubling-remainder` for the closed
complement).

Reproduce the frozen peak decomposition with:

```bash
./build/navier_stokes_lab local-sld-block \
  --state proof/l4/states/local-sld-trajectory/maximum-T050-K3-from-K2/K3.tsv \
  --doubling-family --threads 12 --evolve-steps 596 \
  --nu 0.1 --dt 0.0005 \
  --certificate proof/l4/analysis/shifted-local-density/frozen-maximum-sld-T050-K3-from-K2-peak-blocks.json
```
