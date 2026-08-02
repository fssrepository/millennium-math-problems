# Dyadic projective-height Schur reduction

This note gives an exact algebraic reduction of the open projective part of
RQ-11 to two cutoff-uniform dyadic estimates. It does not prove either
uniform estimate and therefore does not prove RQ-11, the local SLD lemma, or
the Clay problem.

## Exact height matrix

Partition the primitive feasible squared-length shapes by

```text
Gamma_0={c=1},
Gamma_j={2^(j-1)<c<=2^j},  j>=1,
```

where `c` is the largest primitive squared length. Write `B_j` for the sum of
all projective advection operators in `Gamma_j`, and set

```text
b_j=B_j(u,u),  c_j=B_j(u,Au),
S_j=<Au,b_j>,  T_j=<A b_j,Au>.
```

The diagonal quartet `J_jj` and the symmetric off-diagonal quartet `J_ij`
are the five-component polarization of the closed bracket:

```text
J_jj = -<b_j,A b_j> + <b_j,c_j> - <Au,B_j(b_j,u)>
       + S_j^2/(2Z) + 3 S_j T_j/(2P),

J_ij = -<b_i,A b_j>-<b_j,A b_i>
       +<b_i,c_j>+<b_j,c_i>
       -<Au,B_i(b_j,u)>-<Au,B_j(b_i,u)>
       +S_i S_j/Z+3(S_i T_j+S_j T_i)/(2P).          (DHS-1)
```

`LocalSldProjectiveHeightMatrix` evaluates every term in DHS-1 directly. Its
sum reconstructs the selected full bracket below `1e-17` relative error on
all saved stress states.

## Absolute-component Schur lemma

Let

```text
q=|S_full|/(Z^2 P^2).
```

Define `e_j` as `q` times the sum of the absolute values of the five
components of `J_jj`, and define `e_ij` analogously from the five components
of `J_ij`. For every pair with nonzero diagonal weights set

```text
r_jj=1,
r_ij=e_ij/(2 sqrt(e_i e_j)),  i!=j,
R=max_i sum_j r_ij.                                  (DHS-2)
```

Then the arithmetic--geometric mean inequality gives the exact finite-matrix
bound

```text
sum_j e_j + sum_(i<j) e_ij
 <= sum_j e_j + sum_(i<j) r_ij(e_i+e_j)
 <= R sum_j e_j.                                    (DHS-3)
```

The left side dominates the absolute value of the complete projective
power-one block. Therefore RQ-11 follows if both

```text
sup_(N,u) R_N(u) < infinity,                         (DHS-4)
sup_(N,u) sum_j e_j,N(u) < infinity                  (DHS-5)
```

hold on the normalized state class. A sufficient form of DHS-5 is a dyadic
decay `e_j<=C 2^(-alpha j)` for any `alpha>0`. DHS-3 is proved algebraically;
DHS-4 and DHS-5 are the remaining analytic statements.

The component envelope is essential. Normalizing off-diagonal entries by the
absolute value of the signed diagonal `J_jj` is invalid because internal
cancellation can make `J_jj` almost zero. The finite ratio then exceeds
`10^7` on a saved K8 state. The five-component envelope removes that
artificial small denominator.

## Current falsification results

On the H=64 K8 open-power winner, the exact matrix gives

```text
total component envelope      0.0038520661
diagonal component envelope   0.0021599994
maximum Schur row sum          2.9859137
Schur upper bound              0.0064495718
actual / Schur upper bound     0.5972592.
```

The maximum normalized off-diagonal envelope is `0.5372`. By shell gap its
observed maxima decrease as

```text
gap       1       2       3       4       5       6       7       8
ratio   .537    .426    .283    .246    .158    .0529   .00789  .000125.
```

The same maximum row sum stays between `2.89` and `3.85` on the K8
open-power, height-power, and height-stretching adversary winners tested so
far. This is finite evidence for DHS-4, not a uniform theorem.

An exact-gradient objective maximizes the complete signed diagonal block

```text
|J_(H,2H] S_full|^2/(Z^4P^4).
```

On separately optimized K8 states its square roots at H=8,16,32,64 are
approximately

```text
5.47e-4, 4.73e-4, 3.93e-4, 2.45e-4.
```

The fitted signed-diagonal height slope is `-0.374`. More importantly, the
absolute five-component diagonal envelope has fitted slope `-0.435`, so the
observed decay is not caused only by cancellation inside `J_jj`. These are
optimized finite lower branches, not the upper bound DHS-5.

## Reproduction

```bash
./build/navier_stokes_lab local-sld-projective-height-matrix \
  --state proof/l4/states/local-projective-open-power/height-64-K8-warm/K8.tsv \
  --threads 12 --exclude-triple-family \
  --certificate proof/l4/analysis/shifted-local-density/remainder-quartet/K8-open-power-height64-winner-height-matrix.json

./build/navier_stokes_lab local-sld-projective-height-transfer \
  --height 8 --state proof/l4/states/local-projective-height-power/H8-K8-warm/K8.tsv \
  --height 16 --state proof/l4/states/local-projective-height-power/H16-K8-warm/K8.tsv \
  --height 32 --state proof/l4/states/local-projective-height-power/H32-K8-warm/K8.tsv \
  --height 64 --state proof/l4/states/local-projective-height-power/H64-K8-warm/K8.tsv \
  --threads 12 --exclude-triple-family \
  --certificate proof/l4/analysis/shifted-local-density/remainder-quartet/K8-optimized-height-power-transfer.json
```

Artifacts:

- [`../../analysis/shifted-local-density/remainder-quartet/K8-open-power-height64-winner-height-matrix.json`](../../analysis/shifted-local-density/remainder-quartet/K8-open-power-height64-winner-height-matrix.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K8-optimized-height-power-transfer.json`](../../analysis/shifted-local-density/remainder-quartet/K8-optimized-height-power-transfer.json)

