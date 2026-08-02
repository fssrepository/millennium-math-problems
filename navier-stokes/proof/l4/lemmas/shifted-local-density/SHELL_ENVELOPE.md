# Local quartic shell envelope

This note proves the cutoff-independent estimate used by
`LocalQuarticShellEnvelope`. It is a partial lemma for SLD-1P, not a proof of
SLD-1P itself.

Let `A=-Delta`, let `B_L=B_local(u,u)`, and define dyadic shells

```text
Lambda_j = {k: R_j <= |k| < 2 R_j},   R_j=2^j,
E_j      = sum_{k in Lambda_j} |u_k|^2,
E_j^near = E_{j-1}+E_j+E_{j+1}.
```

Missing negative shells are interpreted as zero. The local mask requires the
largest and smallest frequencies in a triad to differ by at most two.

## Per-shell estimate

Fix an output `k` in `Lambda_j`. Every local input satisfies

```text
R_j/2 <= |p|,|q| < 4 R_j.
```

The Leray projector has operator norm one. There are at most
`(8R_j+1)^3 <= 729 R_j^3` possible lattice inputs, and `|q|^2 < 16R_j^2`.
Cauchy--Schwarz in the convolution sum therefore gives

```text
|B_L(k)|^2
 <= 11664 R_j^5 sum_{p+q=k, local} |u_p|^2 |u_q|^2.
```

After summing over `k`, each ordered pair `(p,q)` has one target and both
inputs lie in the three-shell neighborhood. Since `|k|^2<4R_j^2`,

```text
||(A^(1/2) B_L)_j||_2^2
 <= 46656 R_j^7 (E_j^near)^2.                         (LQE-1)
```

The constant is deliberately elementary rather than optimized. It depends on
neither the Galerkin cutoff nor the state.

## Global summation

Use the nonnegative factorization

```text
R_j^7 (E_j^near)^2
 = (R_j^3 E_j^near)(R_j^4 E_j^near).
```

Three-shell overlap gives the exact coarse constants

```text
sum_j R_j^3 E_j^near <= (2^-3+1+2^3) sum_j R_j^3 E_j
                      = (73/8) sum_j R_j^3 E_j,

sum_j R_j^4 E_j^near <= (2^-4+1+2^4) sum_j R_j^4 E_j
                      = (273/16) sum_j R_j^4 E_j.
```

Because `R_j<=|k|` on `Lambda_j`, Cauchy--Schwarz yields

```text
sum_j R_j^3 E_j
 <= (sum_j R_j^2 E_j)^(1/2)
    (sum_j R_j^4 E_j)^(1/2)
 <= sqrt(Z P),

sum_j R_j^4 E_j <= P.
```

For nonnegative sequences, the sum of pointwise products is no larger than
the product of the sums. Summing LQE-1 consequently proves

```text
||A^(1/2) B_L||_2^2
 <= 46656 (73/8)(273/16) sqrt(Z P) P
 = 7264120.5 sqrt(Z) P^(3/2).                         (LQE-2)
```

LQE-2 is uniform in the Fourier cutoff and uses only the current enstrophy and
palinstrophy. In particular, it does not assume a future `H^s`, `s>2`, bound.

## Machine certificate

For every reported shell, the C++ certificate records the actual left side of
LQE-1, its explicit right side, all three geometry ratios, the two overlap
moments, and both stages of the global summation. On the optimized K6 SLD
state:

```text
max per-shell actual/bound        1.2801803470e-5
global actual                     0.603413126167
global shell-product bound        6.01892325668e6
global Z/P bound                  1.04930748081e7
global actual/ZP-bound ratio      5.7505844302e-8
```

The large slack is harmless for cutoff independence, but it means this bound
should not be used to predict the sharp numerical SLD constant.

## Remaining obligation

LQE-2 closes the isolated outer-slot norm without a high-Sobolev assumption.
SLD-1P still contains the advecting variation and nonlinear palinstrophy
variation. Their signed combination with the outer and advected slots must be
estimated before the shifted logarithmic inequality follows.
