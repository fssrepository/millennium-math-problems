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

The saturating state is not opaque. For the normalized cyclic axis shear

```text
u_hat(e1)=e3/sqrt(3),
u_hat(e2)=e1/sqrt(3),
u_hat(e3)=e2/sqrt(3),
u_hat(-k)=conj(u_hat(k)),
```

the engine certifies

```text
E=Z=P=1,  S=0,  K+G=-1/3,  c=-1/3
```

with maximum identity error `1.08e-19`. Padding this same state and running
12 independent starts finds `|c|=1/3` unchanged from K1 through K8, with no
larger finite-dimensional example. This motivates the sharp static candidate

```text
|K+G| E^(1/4) <= (1/3) Z^(7/4) P.                   (LQC-7)
```

LQC-7 is not proved. Even if proved, it would use the current state's
normalization; LQC-3 along a trajectory uses `k0,B0` frozen at time zero.
That distinction prevents the sharp static conjecture from being presented
as the completed local trajectory lemma.

The signed variant does not have a universal favorable sign either. Exact
reverse mode for

```text
c=(K+G)E^(1/4)/(Z^(7/4)P)
```

agrees with central differences to `1.57e-11`. A 12-start search finds the
positive value `c=0.0740187069851`, stable from K2 through K4. The corresponding
stretching is nearly zero, so this rejects the shortcut `K+G<=0` without
creating a large direct SLD source. The same test also finds positive brackets
inside the doubling family and its complement. This failed route is recorded
as F012; subsequent estimates must preserve the joint factorization LQC-6.

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

The instantaneous LQC-5 oracle uses the evaluated state as its datum. A
separate `LocalSldTrajectoryAdjoint` now evaluates the terminal state with
`k0=sqrt(Z(0)/E(0))` and `B0=E(0)P(0)` frozen from the actual initial state.
Its reverse pass includes both the RK4 VJP and the direct dependence of the
two frozen parameters on the optimized initial datum. Central differences
give `5.62e-12`; at zero steps it reconstructs the instantaneous gradient to
`1.52e-19`.

At `nu=0.1`, 12-start K3 horizon continuation gives

```text
T          terminal frozen R_local       dt-halving error
0.01       7.98918e-4                    3.64e-18
0.05       8.12240e-4                    1.39e-17
0.10       8.26646e-4                    9.14e-17
0.20       8.46863e-4                    4.21e-16
```

The observed increase proves that the static extremizer cannot replace a
frozen-data trajectory estimate.
The nonsmooth maximum-on-trajectory objective over `0<=t<=0.5`, warm-started
from the K2 winner, selects `t=0.298` on both the coarse and refined grids and
gives the current K3 lower bound `8.53498799310e-4`. The refined value differs
by `8.16e-16` relatively. This supersedes the earlier `8.48675785e-4`
secondary branch. The projected initial gradient is still `4.04e-4`, so the
new value is not claimed as a converged global maximum. At the refined peak,
the exact frozen-normalization decomposition is `8.24077287430e-4` from the
doubling family, `1.51420260355e-5` from the closed remainder, and
`1.42794858448e-5` from the mixed term. They reconstruct the full value to
`6.20e-20`; the family retains `96.55%` of the absolute three-block sum.
The initial state has `99.969%` of its energy in the first hard shell and the
`(1,1,2)` signature carries `99.627%` of coherent local transfer.

## Structured response hierarchy

The low-shell pattern is reproducible without a full-state optimizer. Starting
from the cyclic axis shear, define orthonormal response states recursively by

```text
raw b_n = sum_{i+j=n-1} B(b_i,b_j),
b_n     = normalized Gram--Schmidt(raw b_n; b_0,...,b_{n-1}).
```

The scalar chain misses one polarization in the `(2,1,1)` orbit and two
oriented `(3,1,0)` cyclic orbits. `LocalSldCyclicOrbitBasis` supplies explicit
divergence-free representatives for those three directions. On the current
K3 winner, orders `0..15` plus the three orbit directions give

```text
captured state energy              0.999998627866
state-space residual norm          0.00117138141911
projected frozen trajectory peak   8.53172846573e-4
unrestricted frozen peak           8.53498799310e-4
captured objective fraction        0.999618098189
peak time for both states           0.298
```

The projected trajectory was evaluated independently by the direct and FFT
RK4 implementations; their refined values agree to all reported digits. This
is strong evidence for a compact extremal response structure. It is not a
proof that every cutoff or every smooth datum stays in this structure. The
analytic target is now a uniform summability bound for the response
coefficients together with a controlled complement estimate.
The next analytic task is to bound this evolving joint quotient uniformly,
using the exact block split below. No finite horizon scan establishes that
bound.

[SIGNATURE_BLOCK.md](SIGNATURE_BLOCK.md) now makes this routing exact. It
splits the local operator into the `(m,m,2m)` family and its complement and
independently reconstructs the closed-family, closed-remainder, and mixed
brackets. On the K8 winner the family contributes `7.70182e-4` of the full
`7.95365e-4` quotient; remainder and mixed contributions are both about
`1.3e-5`.

Exact gradients with the common full stretching factor now maximize the three
contributions independently. The observed maxima are `7.72930e-4` for the
doubling family through K6, `2.20683e-4` for the remainder through K6, and
`1.34936e-4` for the mixed block through K6. Because these winners differ,
the values are not an additive full-state bound. They prove that a successful
argument cannot discard the latter two blocks solely from their small size on
the original K8 winner.

Cutoff stability of `C_state` is a falsification check only. The proof task is
to derive LQC-3, or a weaker signed version sufficient for LQC-4, uniformly
along every Galerkin trajectory issued from one fixed smooth datum.
