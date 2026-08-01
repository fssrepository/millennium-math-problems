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

## Direct signed target

The absolute LQC-3 screen is stronger than the local polynomial lemma. The
engine therefore also maximizes the actual signed quotient at the initial
time:

```text
R_local(u)=4 S^3 ZP(K+G)
 / [k0(S^4 Z^2P+B0 Z^3P^4)],

k0=(Z/E)^(1/2),  B0=EP.                              (LQC-5)
```

Unlike the absolute `C_state`, LQC-5 retains the `S^3` factor and the complete
positive denominator. Its manually derived reverse-mode gradient traverses
every nested local bilinear-advection node. The central-difference regression
has relative error `6.77e-12`.

At `t=0`, LQC-5 has the exact dimensionless factorization

```text
c = (K+G) E^(1/4)/(Z^(7/4)P),
x = S/(E^(1/4)Z^(1/4)P),

R_local = c 4x^3/(1+x^4).                            (LQC-6)
```

Thus the direct target measures a joint `c,x` envelope. Maximizing `|c|`
alone can select a state where the shape factor vanishes. The engine evaluates
both sides of LQC-6 independently; their relative discrepancy is included in
the self-test.

The distinction matters. An 80-step K1 search makes the stronger absolute
ratio converge to exactly `1/3`, but that state has
`S=-1.23e-11`; consequently its real polynomial numerator is negligible. It
is a legitimate sharp test of LQC-3, not the worst case for SLD-1P-L.

Twelve-start projected L-BFGS optimization of LQC-5 gives

```text
K       R_local                 continuation gain
1       7.604090790190e-4       -
2       7.953253043740e-4       2.410e-5
3       7.953641813777e-4       3.787e-8
4       7.953648543733e-4       3.440e-10
5       7.953649287667e-4       7.439e-11
6       7.953649534336e-4       2.467e-11
7       7.953649580969e-4       4.663e-12
8       7.953649600919e-4       1.995e-12
```

At K8 the numerator and denominator are `0.00221536043169` and
`2.78533823195`; `S=-0.131915633194`, `E=1`, `Z=1.05843710262`, and
`P=1.17949800570`. The K8-to-K7 projection residual is `1.99e-6`. This is
strong cutoff-convergence evidence for one optimized branch, not a uniform
inequality.

## Symmetry-reduced pattern

The K8 winner puts `99.979%` of its energy in the first hard shell. Six nearly
equal axis modes and twelve nearly equal face-diagonal modes dominate, and the
local squared-frequency signature `(1,1,2)` supplies `98.8%` of the signed
enstrophy transfer.

`LocalSldCyclicAnsatz` turns that observation into an explicit two-basis
family: a normalized cyclic axis shear plus its normalized quadratic local
advection response. A one-dimensional C++ search gives

```text
axis energy fraction       0.943691414498
response energy fraction   0.056308585502
R_local                    0.000775010869995
```

This captures `97.44%` of the full K8 adversarial value. The nonzero full
projected gradient (`0.00280`) measures the missing higher-response
correction, rather than hiding it.

LQC-5 as implemented here is an instantaneous `t=0` oracle: the evaluated
state supplies `k0` and `B0`. The proof still needs a bound at later times
with those two quantities frozen at their actual initial values. The next
analytic task is to prove the signed `(1,1,2)` block estimate and bound its
orthogonal shell remainder with a cutoff-independent summation, then extend
the estimate along a Galerkin trajectory. No finite table establishes those
steps.

Cutoff stability of `C_state` is a falsification check only. The proof task is
to derive LQC-3, or a weaker signed version sufficient for LQC-4, uniformly
along every Galerkin trajectory issued from one fixed smooth datum.
