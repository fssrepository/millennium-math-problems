# Exact Young closure target

The two-entry identity isolates one explicit sufficient estimate for the local
part of SLD-1P. With `K` and `G` from `REDUCED_QUARTET.md`, it is

```text
|K+G| <= C k0 B0^(1/4) Z^(5/4) P^(3/4).             (LQC-3)
```

This file proves the algebraic implication from LQC-3 to SLD-1P-L. It does
not claim that a cutoff-independent `C` has been proved.

Define the two positive denominator terms

```text
X=S^4 Z^2 P,
Y=B0 Z^3 P^4.
```

Their weighted geometric mean is

```text
X^(3/4)Y^(1/4)
 = |S|^3 B0^(1/4) Z^(9/4) P^(7/4).
```

Multiplying LQC-3 by `4|S|^3 ZP` therefore gives

```text
4|S|^3 ZP |K+G|
 <= 4 C k0 X^(3/4)Y^(1/4)
 <= C k0 (3X+Y)
 <= 3 C k0 (X+Y).                                   (LQC-4)
```

The middle inequality is weighted AM--GM. Consequently LQC-3 implies the
local polynomial lemma with `A_local=3C`.

Relative to the elementary quartic scale `Z^(1/2)P^(3/2)`, LQC-3 asks for a
`(Z/P)^(3/4)` frequency depletion, compensated by the fixed-data factor
`k0 B0^(1/4)`. This is the precise analytical gain still missing; it is no
longer described only as an unspecified cancellation.

`LocalQuarticClosureTarget` certifies LQC-4 and measures the finite-state ratio

```text
C_state = |K+G|/(k0 B0^(1/4) Z^(5/4)P^(3/4)).
```

For the optimized K6 state:

```text
K+G                               -0.179749427015
C_state                            0.127982322632
certified local constant 3C        0.383946967896
geometric identity error           8.00e-20
Young left/right ratio             0.005367658904
```

The independently optimized cutoff winners give

```text
K             3          4          5          6
C_state    0.145045   0.130869   0.129389   0.127982
```

Zero-padding the identical K6 state to K7 and K8 changes the ratio from
`0.127982322632` to `0.127982322669` and `0.127982322689`. This branch is
cutoff stable, but the table remains a falsification result rather than the
uniform estimate required by LQC-3.

Cutoff stability of `C_state` is a falsification check only. The proof task is
to derive LQC-3, or a weaker signed version sufficient for LQC-4, uniformly
along every Galerkin trajectory issued from one fixed smooth datum.
