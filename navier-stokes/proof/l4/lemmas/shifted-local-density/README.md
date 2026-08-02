# Shifted local critical-density candidate

Companion response-space proof targets:

- [`RESPONSE_DIAGONAL.md`](RESPONSE_DIAGONAL.md) defines the boundary-free
  cutoff diagonal and its weighted coefficient majorant.
- [`RESPONSE_TENSOR.md`](RESPONSE_TENSOR.md) records the exact interaction
  tensor and the explicit two-radius bilinear candidate that survived every
  ordered response pair through K12.
- [`DOUBLING_QUARTET.md`](DOUBLING_QUARTET.md) proves the half-derivative
  one-shell gain for the dominant orthogonal family and identifies the exact
  signed cross-shell obstruction still blocking its global sum.

This note records a scale-compatible replacement for the rejected pointwise
monotonicity and unshifted multiplicative-growth routes.

Define

```text
C_N(t) = |V_N,local(t)|^4/(Z_N(t) P_N(t)^3),
B_0    = E_N(0) P_N(0),
k_0    = sqrt(Z_N(0)/E_N(0)).
```

The exact candidate is

```text
L4.2-S:
d/dt log(C_N(t)+B_0) <= A(u_0,nu,T) k_0 Z_N(t)              (SLD-1)
```

with `A` independent of the Galerkin cutoff. SLD-1 is open.

## Exact compatibility checks

Under velocity-amplitude scaling, both `C_N` and `B_0` have degree four.
Under Navier--Stokes spatial scaling their exponent is two:

```text
C_N: 4*3 - 1 - 3*3 = 2,
B_0: -1 + 3 = 2.
```

The logarithmic time derivative has scaling exponent two. Since
`k_0` has exponent one and `Z_N` has exponent one, the right-hand side of
SLD-1 also has exponent two. `ShiftedCriticalDensityLemma` verifies these
relations with exact rational arithmetic in every self-test.

## Conditional closure

If SLD-1 is proved, the Galerkin energy identity gives

```text
integral_0^T Z_N(t) dt <= E_N(0)/(2 nu),

C_N(t)+B_0
  <= (C_N(0)+B_0) exp(A k_0 E_N(0)/(2 nu)).
```

Consequently `sup_N integral_0^T C_N(t) dt` is finite for each fixed smooth
initial datum, provided every constant in SLD-1 is cutoff independent. This
would close the remaining local L4.2 block when combined with the proved
moving far-tail estimate and a compatible transition-band estimate.

## Numerical screen

The exact discrete adjoint optimizes

```text
log((C_N(T)+B)/(C_N(0)+B))
```

for a configured finite `B`. At `E(0)=1`, `H4^2<=100`, `nu=0.1`, and
`T=0.001`, the K6 shift sweep gives

```text
B             1e-5       1e-4       1e-3
gain         0.011141    0.004508    0.000938
gain/(T k0 Z) 10.657       4.135       0.826
C_N(0)       1.87e-5     9.74e-5     1.89e-4
Delta C_N    3.21e-7     8.92e-7     1.12e-6
```

Time-step errors are below `2.7e-14`. The dependence on `B` is substantial,
so a free numerical shift cannot be promoted to a lemma. `B_0=E(0)P(0)` is
the current analytical choice because it has the exact required homogeneity
and depends only on the fixed smooth initial datum. The next proof task is to
differentiate `C_N` analytically and decide whether its local triad terms can
be bounded by the right-hand side of SLD-1 without a future high-Sobolev norm.

## Exact state-dependent shift adversary

The optimizer also differentiates `B_0=E(0)P(0)` as part of the initial state;
it is not frozen as a numerical parameter. The extra gradient term is

```text
dB_0 = P(0) dE(0) + E(0) dP(0).
```

The complete discrete gradient agrees with centered differences to
`2.06e-10`. Twelve starts and eight L-BFGS iterations give

```text
K                                  3          4          5          6
log((C(T)+B0)/(C(0)+B0))      7.918e-7   8.310e-7   8.558e-7   8.632e-7
gain/(T k0 Z(0))              7.240e-4   7.533e-4   7.751e-4   7.810e-4
```

Because different restarts win these four searches, the K6 winner was then
replayed independently at every lower cutoff. Finally, the identical K6 state
was zero-padded to K7 and K8 without optimization:

```text
K                                  6             7             8
shifted log gain              8.632216851e-7 8.632152689e-7 8.632152671e-7
gain/(T k0 Z(0))              7.810352793e-4 7.810294739e-4 7.810294723e-4
relative change                    -        -7.43e-6      -2.14e-9
```

This finite branch is time-step stable and cutoff converged. It does not prove
SLD-1 or bound its global constant, but it does not falsify the candidate.

## Instantaneous derivative oracle

`shifted-density` evaluates the left side of SLD-1 directly, without an RK4
horizon:

```text
dC_local/dt = <gradient C_local, Navier--Stokes RHS>,
R_SLD       = (d/dt log(C_local+E0P0))/(k0 Z(0)).
```

The original scalar oracle evaluates the K6 winner in `0.22 s`; the complete
12-thread derivative/role ledger takes `0.48 s`. Direct-triad and FFT backends
give identical main diagnostic values, while the expanded budget components
agree within `1.1e-19`. The separately optimized winners give

```text
K                  3             4             5             6
R_SLD       7.239567e-4  7.534130e-4  7.754572e-4  7.815121e-4
```

At K6 the instantaneous rate differs from the `T=0.001` average by only
`6.10e-4` relatively. This command is now the fast falsification path for
SLD-1 states; trajectory integration is reserved for candidates that survive
the instantaneous screen.

Reproduce with

```bash
./build/navier_stokes_lab shifted-density \
  --state proof/l4/states/local-critical-ep-log-gain-multistart/dynamic/K6.tsv \
  --nu 0.1 --threads 12 --backend fft \
  --certificate /tmp/K6-shifted-density.json
```

Artifacts:

- [`../../adversary/shifted-local-density/local-critical-shifted-log-gain1e-5-h4-cap100-K6-multistart6-lbfgs8.json`](../../adversary/shifted-local-density/local-critical-shifted-log-gain1e-5-h4-cap100-K6-multistart6-lbfgs8.json)
- [`../../adversary/shifted-local-density/local-critical-shifted-log-gain1e-4-h4-cap100-K3-K6-multistart12-lbfgs8.json`](../../adversary/shifted-local-density/local-critical-shifted-log-gain1e-4-h4-cap100-K3-K6-multistart12-lbfgs8.json)
- [`../../adversary/shifted-local-density/local-critical-shifted-log-gain1e-3-h4-cap100-K6-multistart6-lbfgs8.json`](../../adversary/shifted-local-density/local-critical-shifted-log-gain1e-3-h4-cap100-K6-multistart6-lbfgs8.json)
- [`../../adversary/shifted-local-density/local-critical-ep-log-gain-h4-cap100-K3-K6-multistart12-lbfgs8.json`](../../adversary/shifted-local-density/local-critical-ep-log-gain-h4-cap100-K3-K6-multistart12-lbfgs8.json)
- [`../../adversary/shifted-local-density/local-critical-ep-log-gain-K6-projections-K3-K6.json`](../../adversary/shifted-local-density/local-critical-ep-log-gain-K6-projections-K3-K6.json)
- [`../../adversary/shifted-local-density/local-critical-ep-log-gain-K6-lifts-K6-K8.json`](../../adversary/shifted-local-density/local-critical-ep-log-gain-K6-lifts-K6-K8.json)
- [`../../analysis/local-critical-increase/K6-instantaneous-E0P0-shifted-density.json`](../../analysis/local-critical-increase/K6-instantaneous-E0P0-shifted-density.json)

## Exact derivative ledger

Write `S=V_N,local` with its sign, so `C=S^4/(Z P^3)`. The new
`LocalCriticalDerivativeLedger` evaluates the exact chain rule

```text
C' = 4 S^3 S'/(Z P^3)
     - S^4 Z'/(Z^2 P^3)
     - 3 S^4 P'/(Z P^4).
```

Every term is independently split into the Euler and viscous parts of the
Galerkin right-hand side. The reconstructed derivative agrees with the
independent `gradient C dot RHS` oracle to `1.18e-18` relatively on the K6
winner. The ledger also certifies

```text
Z'_NL = -2 S_global,
Z'_nu = -2 nu P,
P'_nu = -2 nu H3,
H3    = sum_k |k|^6 |u_k|^2.
```

For `S=<A u,B_local(u,u)>`, `StretchingDerivativeLedger` keeps the three
Frechet slots separate:

```text
S'[h] = <A h,B_local(u,u)>
      + <A u,B_local(h,u)>
      + <A u,B_local(u,h)>.
```

At K6, for the nonlinear direction `h=-B_all(u,u)`, these terms are

```text
outer state       -0.609617445980
advecting slot    +0.016052994246
advected slot     +0.328901428531
total             -0.264663023203
```

Their sum agrees with the VJP derivative to `9.22e-19` relatively. For the
viscous direction, the advected slot vanishes to `1.03e-20` relative to the
other two slots. This is the transport cancellation
`<A u,B_local(u,A u)>=0`, preserved because the local triad mask is symmetric.
`LocalQuarticIdentityLedger` further certifies the exact signed square

```text
<A(-B_local),B_local> = -||A^(1/2) B_local||_2^2.
```

Both sides equal `-0.603413126167` at K6 with zero serialized residual. Thus
the dominant outer term has a rigorous sign and norm representation; only its
combination with the other slots and the `P'` contribution remains open.

The dimensionless K6 SLD-1 budget is

```text
source                         contribution /(k0 Z)
S' chain term                   +7.9252418586e-4
Z' chain term                   -2.3985915123e-6
P' chain term                   -8.6135345824e-6
local nonlinear RHS             +6.8716890362e-4
nonlocal nonlinear RHS          +6.4812110181e-5
viscous RHS                      +2.9531045965e-5
total                            +7.8151205976e-4
```

Thus `91.38%` of the nonlinear source on this extremizer is the pure local
quartic block. The nonlocal feedback is measurable but secondary. Viscosity
does not have a fixed favorable sign for the quotient: at K6 its net
contribution is positive even though it dissipates both `Z` and `P`.

The pure local normalized block retains substantial internal cancellation:

```text
outer negative square             +2.3068042175e-3
advecting slot                    -6.1297970097e-5
advected slot                     -1.2458880655e-3
local enstrophy term              -3.5140539404e-5
local palinstrophy term           -2.7730873889e-4
reconstructed local source        +6.8716890362e-4
```

The last four terms cancel `70.21%` of the outer source. The typed budget
reconstructs the local nonlinear value to `6.93e-19` relatively. This is why
the proof target is the combined quartet expression, not five independent
absolute-value estimates.

Zero-padding the identical K6 state gives total instantaneous rates
`7.8151205976e-4`, `7.8151207993e-4`, and `7.8151208069e-4` at K6, K7, and K8.
Only the nonlinear stretching channel changes. The K7-to-K8 increment is
`7.54e-13` absolutely, so the role decomposition itself is cutoff converged
on this branch.

## Shell-resolved local quartet

`LocalQuarticShellLedger` groups the complete local-local contribution by
target mode, exact Laplace eigenshell, hard cutoff shell, and dyadic radial
annulus. On the optimized K6 state, cancellation within a target mode is
`47.13%`; aggregation within one dyadic annulus raises this to `70.66%`, while
cancellation between dyadic annuli is only `3.21e-7`. The normalized annular
totals are

```text
j=0   +6.2714763447e-4
j=1   +6.0016171958e-5
j=2   +5.2074922362e-9
j=3   -1.1029824132e-10.
```

This localization is state dependent. The critical-integral and
local-signature control states instead cancel about `99%` between their first
two dyadic annuli. Therefore no proof step assumes that the observed K6
within-annulus cancellation is universal; signed shell totals or their
positive parts must be controlled.

The isolated negative-square norm now has a conventional cutoff-independent
shell estimate. For `R_j=2^j` and three-shell neighborhood energy
`E_j^near`, `LocalQuarticShellEnvelope` certifies

```text
||(A^(1/2)B_local)_j||_2^2
 <= 46656 R_j^7 (E_j^near)^2,

||A^(1/2)B_local||_2^2
 <= 7264120.5 sqrt(Z) P^(3/2).                       (LQE-2)
```

The second line follows by factoring `R^7 E^2=(R^3E)(R^4E)`, bounded
three-shell overlap, and `H^(3/2)` interpolation between `Z` and `P`. It uses
no future norm above `H^2`. The K6 global actual-to-bound ratio is
`5.75e-8`; the size of that slack is irrelevant to cutoff independence.
See [SHELL_ENVELOPE.md](SHELL_ENVELOPE.md) for the human-checkable proof.

`LocalQuarticCommutator` also combines the outer and advected slots exactly:

```text
-<A B_L,B_L>-<A u,B_L(u,B_L)>
  = -<B_L,A B_L-B_L(u,A u)>.
```

The Fourier symbol is `|k|^2-|q|^2` and obeys
`abs(|k|^2-|q|^2)<=|p|(|k|+|q|)`. At K6 the two raw values
`-0.603413126167` and `+0.325899010744` combine to
`-0.277514115423`; the independently evaluated identity has relative error
`1.95e-19`. See [COMMUTATOR.md](COMMUTATOR.md).

There is a second exact grouping that incorporates both denominator
derivatives. With `h=-B_L`, the outer, enstrophy, and palinstrophy entries are

```text
4S^3/(ZP^3) <A h,
    B_L-Su/(2Z)-3S A u/(2P)>.
```

`LocalQuarticProjectedResidual` reconstructs this pairing both by direct
evaluation and by the completed square

```text
-||A^(1/2)(B_L-3S A u/(4P))||_2^2
  + S^2/(2Z) + 9S^2 H3/(16P^2).
```

At K6 the normalized projected, advecting, and advected entries are
`+1.99435e-3`, `-6.12980e-5`, and `-1.24589e-3`; their sum is the complete
`+6.87169e-4` local source. This reduces five terms to three with
`3.85e-18` relative error. The prefactor contains signed `S^3`, so the square
does not have one universal favorable sign. See
[PROJECTED_RESIDUAL.md](PROJECTED_RESIDUAL.md).

Finally, combining the commutator and projected-residual identities reduces
the full local-local polynomial numerator to

```text
K = -<B_L,A B_L-B_L(u,A u)>
    +S^2/(2Z)+3S<A B_L,A u>/(2P),
G = <A u,B_L(-B_L,u)>,

N_local = 4S^3 ZP (K+G).                              (LQR-4)
```

`LocalQuarticReducedLedger` checks this directly against the original
five-entry and denominator-free calculations. At K6, `K` and `G` contribute
`+7.48467e-4` and `-6.12980e-5` after normalization. The polynomial local
numerator is `0.00264627230688`, with `3.20e-19` relative reconstruction
error. [REDUCED_QUARTET.md](REDUCED_QUARTET.md) states the resulting two-entry
local lemma `SLD-1P-L`.

The two-entry formula yields a concrete sufficient target:

```text
|K+G| <= C k0 B0^(1/4) Z^(5/4)P^(3/4).               (LQC-3)
```

Weighted AM--GM then gives the local SLD-1P bound with `A_local=3C` exactly.
The required improvement over the elementary quartic scale is a
`(Z/P)^(3/4)` depletion. `LocalQuarticClosureTarget` certifies the algebra and
reports `C_state=0.127982322632` at K6; this finite ratio is not a proof of a
uniform constant. The K3--K6 optimized values decrease from `0.145045` to
`0.127982`, and the identical K6 state changes by only `5.7e-11` after lifting
to K8. See [CLOSURE_TARGET.md](CLOSURE_TARGET.md).

The search now also differentiates the actual signed local SLD-1P-L quotient,
so the stronger absolute target cannot distract the optimizer with states
where `S^3` vanishes. Its K1--K8 continuation converges from
`7.60409e-4` to `7.95364960e-4`; the last continuation gain is `2.0e-12`.
The K8 state is dominated by the local `(1,1,2)` signature. An explicit cyclic
axis-shear plus quadratic-response ansatz recovers `97.44%` of that value.
These results sharpen the next signed block estimate but do not prove it.
The operator-level version is in
[SIGNATURE_BLOCK.md](SIGNATURE_BLOCK.md): it evaluates the dominant closed
block, the closed complement, and every mixed bilinear term independently.
Their K8 reconstruction error is `2.20e-18`.
The corresponding common-normalization gradients pass central differences at
`1.29e-11` or better. Independent scans reach `7.72930e-4`, `2.20683e-4`,
and `1.34936e-4` for the family, remainder, and mixed terms respectively.
The bracket alone is not sign-definite: a replayable search gives positive
normalized `K+G` equal to `0.0740187069851`, while its nearly zero stretching
suppresses the actual SLD contribution. The active lemma must therefore keep
the joint block factorization.

The pure cyclic axis state gives the exact sharp static candidate `|c|=1/3`
with `1.08e-19` identity error, unchanged under padded K1--K8 searches. The
trajectory engine now freezes `k0,B0` at time zero and differentiates through
RK4 plus both parameter dependencies. Its gradient error is `5.62e-12`.
Optimized K3 terminal ratios increase from `7.98918e-4` at `T=0.01` to
`8.46863e-4` at `T=0.20`, so the unresolved lemma is explicitly the uniform
frozen-data trajectory bound rather than the static conjecture alone.
Warm continuation from K2 corrects that older secondary branch. K3 selects
`t=0.298` and reaches `8.53498799310e-4`; six accepted K4 steps improve the
finite lower bound to `8.53527437357e-4` at the same time. The K4 dt-halving
error is `9.00e-16`, hard-shell-four energy is `1.97e-9`, and the final
projected objective gradient is `5.17e-5`. The search is not proved globally
stationary, so this remains a pattern and falsification artifact.
At the refined K4 peak, frozen-data block analysis gives
`8.24057223723e-4` from the doubling family, `1.51705418626e-5` from the
closed remainder, and `1.42996717717e-5` from the mixed block. The relative
reconstruction error is `1.36e-18`. Dominance of the same family therefore
survives the evolution rather than being only a `t=0` artifact.

Nineteen explicit response/orbit directions capture `99.9997617%` of the K4
state energy and `99.9331554%` of its trajectory objective. Fixed-depth
response comparisons beyond order K are polluted by the Galerkin boundary,
so [RESPONSE_DIAGONAL.md](RESPONSE_DIAGONAL.md) retains only orders `0..K`
and states the weighted response-majorant target. The observed diagonal is
evidence, not the required cutoff-uniform bilinear and complement estimate.

## Reduced proof obligation

Multiplying SLD-1 by the positive denominator `Z^2 P^4` removes all quotient
singularities. The remaining pointwise polynomial inequality is

```text
4 S^3 S' Z P - S^4 Z' P - 3 S^4 Z P'
  <= A k0 (S^4 Z^2 P + B0 Z^3 P^4).                 (SLD-1P)
```

This is algebraically equivalent to SLD-1 wherever `Z,P>0`, and it extends
continuously through `S=0`. The active analytical task is now to symmetrize
the local-local quartic contributions to `S'` together with the `P'` term;
bounding the three stretching slots independently would discard the observed
large cancellation and reintroduce a non-closing high-Sobolev norm. The
nonlocal part can be routed through the existing moving far-tail estimate,
leaving the finite transition block explicit.

`ShiftedCriticalDensityBudgetAnalyzer` evaluates SLD-1P through this second,
denominator-free route. At K6 its polynomial numerator is
`0.00300958572245`, its positive denominator is `3.85097796617`, and their
ratio is `7.81512059765e-4`. The relative error against the original shifted
logarithmic route is `4.74e-19`.
