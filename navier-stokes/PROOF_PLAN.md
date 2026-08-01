# Proof roadmap and restart point

Last updated: 2026-08-01

## Exact target

Prove statement **(B)** in Charles Fefferman's official Clay problem
description: for every smooth, periodic, divergence-free initial velocity on
`R³/Z³`, viscosity `ν > 0`, and zero external force, there are periodic
`C∞` velocity and pressure fields solving Navier–Stokes on all of
`R³ × [0,∞)`.

Proving (B) alone is a solution of the Millennium problem. Transferring the
result to `R³` is not required by the official disjunction. The code uses a
`2π` period because it makes Fourier wave numbers integral; this is converted
to unit period by a fixed rescaling.

Numerics are only fast lemma falsifiers and algebra checks. They are not an
accepted terminal step.

## Proof chain

| ID | Statement | Status | Machine task |
|---|---|---|---|
| L0 | Galerkin solutions and the cutoff-uniform energy inequality | known | regression check |
| L1 | Fourier-triad energy cancellation and vortex-stretching identity | known; encoded | detailed triad identity check |
| L2 | Standard monomial `E-Z-P` estimates cannot close the proof | certified in its stated scope | exact scaling/homogeneity certificate |
| L3 | Decompose nonlinear transfer into local and nonlocal frequency flux | encoded; analytical bounds open | exact triad classification and flux ledger |
| L4 | **Key lemma:** a dynamic/geometric flux bound, uniform in Galerkin cutoff and time-integrable in a critical class | open | generate candidates and immediately falsify invalid ones |
| L5 | L4 plus energy inequality implies a cutoff-uniform critical regularity bound | open | symbolic exponent and constant check |
| L6 | Compactness limit, global smoothness, pressure reconstruction | standard after L4/L5, details open | assumption and limit certificate |

## Acceptance conditions for L4

A candidate advances only if all conditions hold:

1. correct amplitude homogeneity and high-frequency Navier–Stokes scaling;
2. no counterexample among exact finite Fourier triads;
3. constant independent of the Galerkin frequency cutoff;
4. right-hand side is time-integrable from previously proved bounds;
5. no hidden assumption of the regularity being proved;
6. a reproducible machine-readable certificate is emitted;
7. after computational screening, every infinite-dimensional estimate has a
   conventional human-checkable proof.

## Work cycle

1. State one lemma candidate as an explicit formula.
2. Run homogeneity and scaling checks.
3. Run detailed Fourier-triad counterexample search.
4. If it survives, derive a proof decomposition and record every dependency.
5. If it fails, append its first concrete obstruction to
   `proof/failed_lemmas.tsv`.
6. Always resume at the lowest-numbered open lemma.

## Current result and next action

`L2` considers

```text
E = ||u||₂²,  Z = ||∇u||₂²,  P = ||Δu||₂²,
|V(u)| <= C E^a Z^b P^c.
```

Exact rational constraints from cubic amplitude homogeneity and 3D scaling
force the best absorbable estimate to leave at least a `Z³` term. This does not
exclude finite-time growth of enstrophy. See
`proof/l2/l2-certificate.json`.

The current restart point is **L3 → L4**: use the encoded local/nonlocal flux
partition to formulate a depletion quantity that controls vortex stretching,
then reject it unless its bound is cutoff-independent and closes L5.

The first explicit L4 candidate is now

```text
D_N(t) = |V_N(t)| / (Z_N(t)^(3/4) P_N(t)^(3/4)),
sup_N integral_0^T D_N(t)^4 Z_N(t)^2 dt < infinity  for every finite T.
```

This is sufficient because Young's inequality gives a Gronwall coefficient
`D_N^4 Z_N^2`. A stronger route using only the energy-level time integral of
`Z_N` would require a quarter-power depletion. The scaling engine verifies the
`1/4` exponent exactly for every absorbable `E^a Z^b P^c` candidate. This is a
requirement, not yet a proof of the candidate bound.

The stronger pointwise proposal `sup_N D_N^4 Z_N < infinity` is rejected. Exact
localized scaling followed by fixed-energy normalization gives
`Q_lambda = lambda^2 Q`. The time-integrated L4-A quantity is different:
`D^4 Z^2` has scaling exponent `+2` and `dt` has exponent `-2`, so its integral
is exactly critical and survives this obstruction.

This rejection concerns a universal bound over all fixed-energy Fourier
states. For one fixed smooth initial datum, the trajectory-restricted target

```text
L4-S: sup_N sup_0<=t<=T Q_N(t) < infinity,  Q_N=D_N^4 Z_N,
```

is a stronger sufficient sublemma that remains open. The exact factorization
`D_N^4 Z_N^2=Q_N Z_N` and the Galerkin energy identity give

```text
integral_0^T D_N^4 Z_N^2 dt <= (sup Q_N) E_N(0)/(2 nu).
```

The engine now certifies this reduction and directly optimizes/measures the
trajectory maximum of `Q_N` on nested analytic projective families.

The dynamic adversary no longer relies on random mutation as its primary
search direction. A hand-written discrete adjoint differentiates the exact RK4
map and the dealiased pseudospectral nonlinearity. The implementation includes
direct JVP/VJP oracles, FFT JVP/VJP kernels, checkpointed reverse propagation,
the exact gradient of `Q=D^4 Z`, and Riemannian projection onto the fixed-energy
sphere. Central-difference and adjoint-duality tests currently agree at roughly
`10^-13` and `10^-19`, respectively.

This stronger search changes the computational evidence, not the logical
status of L4-S. On a short equal-parameter `K=1,2` comparison, three gradient
steps found trajectory maxima about 130 and 275 times larger than three random
mutation attempts. The first FFT-path `K=5` certificate contains 1,330 modes
and completed two accepted gradient steps. These are finite extremizers to
analyze, not a cutoff-independent rejection or proof.

`L3` experiments now split `V_N = V_N^local + V_N^nonlocal`, where a triad is
local when its largest and smallest wave-number magnitudes differ by at most a
factor of two. The next sublemmas are:

```text
L4.1  uniformly control integral |V_nonlocal|^4 / (Z P^3) dt;
L4.2  uniformly control integral |V_local|^4 / (Z P^3) dt.
```

The scalar inequality `|a+b|^4 <= 8(|a|^4+|b|^4)` then recovers L4-A. Current
adversarial paths indicate that L4.2 is the dominant term; this observation is
only routing information, not a proof.

The exact discrete adjoint now differentiates the actual L4-A objective,

```text
J_N(u_0) = integral_0^T |V_N|^4 / (Z_N P_N^3) dt
         = integral_0^T D_N^4 Z_N^2 dt,
```

including every RK4 stage and the trapezoidal quadrature. Its directional
derivative agrees with central differences to `1.12e-12` on the deterministic
self-test.

The first fixed-energy, initial-`H4^2 <= 100` continuation at `nu=0.1` and
`T=0.01` gives

```text
K:       1          2          3          4 ...      8
J_N: 3.658e-6   5.892e-6   6.000e-6   6.006e-6   6.010e-6.
```

A further 24-step `K=8` warm continuation reaches `6.04549e-6`. The refined
state has `H4^2=15.7589`, energy-balance residual `8.3e-12`, active modes only
through shell three at the `1e-8` relative-energy threshold, and cutoff-shell
energy `1.49e-11`. The finite family therefore approaches a smooth low-mode
branch rather than displaying short-time cutoff concentration. This does not
prove L4-A: it covers one normalization, viscosity, short horizon, and one
optimized branch, while the desired constant may depend on the fixed smooth
datum and finite time but must not depend on `N`.

The immediate analytical task is now sharper: determine whether the observed
low-shell branch can be bounded by a frequency-envelope estimate whose
constant depends only on the fixed initial Sobolev norm, energy, viscosity, and
`T`, without inserting an a priori high-Sobolev bound on the solution. The
immediate computational checks are continuation in `T`, decreasing viscosity,
and independent multistarts of the exact L4-A objective. Any growth branch must
then be classified into the existing local/nonlocal L4.1/L4.2 ledger before it
can affect the proof status.

Two warm continuations did not change that short-time picture. At `K=8`,
doubling the horizon to `T=0.02` produced `J_8=1.20543e-5`, approximately twice
the `T=0.01` value, while `Q(T)/Q(0)=0.9852`. Reducing viscosity from `0.1` to
`0.02` at the same horizon produced `J_8=1.21020e-5` and
`Q(T)/Q(0)=0.9879`. Thus neither check exposed accelerating `Q` growth; longer
horizons and independent branches remain necessary.

The optimizer now emits a complete per-iteration trace and supports projected
L-BFGS. On the same `K=8` warm state, 16 L-BFGS iterations produced about
`34.5` times the objective gain of 16 normalized steepest-gradient iterations.
This exposed an optimizer artifact in the earlier branch: projecting the best
`K=8` state to `K=3` and reoptimizing produced a better low-mode state than the
direct high-cutoff search.

Lifting that state successively through `K=3,...,8` gives

```text
K:       3          4          5          6          7          8
J_K: 6.118399e-6 6.119727e-6 6.119747e-6 6.119751e-6 6.119751e-6 6.119752e-6
```

The corresponding squared initial `H4` norms remain between `13.76` and
`13.88`. Top-shell energy decreases from `8.12e-5` at `K=3` to `1.44e-13` at
`K=8`; the maximum adjacent projection residual is `2.22e-4`. Thus this branch
is consistent with a smooth projective limit, and its cutoff gain is already
below `2e-11` after `K=4`. The evidence rejects neither L4-A nor global
regularity. On this total-objective branch the measured nonlocal fourth-power
integral is roughly three orders of magnitude below the local term, but the
separate nonlocal adjoint shows that this is not a universal statement.

The engine now differentiates `J_local` and `J_nonlocal` separately. Their
directional gradients agree with central differences to `4.23e-12` and
`3.25e-12`. At `K=3`, `E(0)=1`, `nu=0.1`, `T=0.01`, and initial
`H4^2 <= 100`, maximizing `J_local` gives `J_local=4.38853e-6` and essentially
zero `J_nonlocal`; maximizing `J_nonlocal` instead gives
`J_nonlocal=2.10177e-7` and drives the initial state close to the H4 boundary.
The latter state places `9.26%` of its energy above shell one, compared with
`0.2645%` for the local extremizer.

A seven-point initial-H4 cap sweep from `25` through `400` initially produced a
log-log slope `0.889697878`, but this value is rejected as an optimizer artifact.
Longer optimization raised the cap-100 objective by about `7.7x`; the other cap
points were not converged to the same standard. The sweep cannot currently
support a cap exponent. The explicit analytical candidate remains

```text
L4.1-H:
sup_N integral_0^T |V_N,nonlocal|^4/(Z_N P_N^3) dt
    <= C(E(0), nu, T, ||u_0||_H4).
```

The constant may depend on the fixed smooth datum but not on cutoff. The proof
must propagate an initial frequency envelope using only non-circular estimates;
assuming a uniform high-Sobolev solution bound would assume the desired
regularity.

Optimizer progress was separated from cutoff growth by repeatedly projecting a
high-cutoff winner down to `K=3`, finishing the low-mode optimization, and then
lifting it one cutoff at a time. The refined nonlocal values are

```text
K:              3             4             5
J_nonlocal: 1.61848084e-6  1.62065350e-6  1.62074608e-6
increment:       -          2.1727e-9      9.2579e-11
top shell E:  6.2453e-4     1.0462e-6      4.5120e-9
```

The K4-to-K5 relative gain is `5.71e-5`, and direct projection checks show that
the earlier apparent K3-to-K6 growth was almost entirely continued low-mode
optimization. This branch is consistent with a smooth cutoff limit, but it is
not a proof or a global maximum. Computationally, the next checks are
independent multistarts and continuation in `T` and `nu` for the separately
optimized partition. Analytically, the next step is a dyadic paraproduct ledger
that isolates the low-high-high term and states exactly which summable
frequency envelope would close its time integral.

The dyadic ledger sharpens that analytical split. It defines gap zero by
`high/low <= 2`, gap one by `2 < high/low <= 4`, and the far tail by
`high/low > 4`. On the converged full-nonlocal K5 branch, signed gap-one
stretching is `-0.371413`, while the entire far tail is `-8.92578e-5`. Thus the
old nonlocal class was dominated by its first transition band, not by strongly
separated scales. The ledger further records whether the low mode is the
advecting, advected, or target factor; the first two roles dominate and the
target role partially cancels them on this branch.

The engine now exposes exact `critical-near-nonlocal-integral` and
`critical-far-nonlocal-integral` objectives. The nonzero far-tail gradient
agrees with central differences to `4.84e-13`. A far-tail adversary validated
with `dt=0.0005` and `dt/2` gives

```text
K:          3             4             5             6
J_far: 2.11194e-11   3.20285e-11   3.21555e-11   3.21556e-11
top E:  1.05881e-3    2.68442e-5    6.81316e-8    1.37352e-10
```

These are separate projections of the same K6 winner, so the saturation is not
an optimizer-progress artifact. The revised analytical ledger is therefore:

```text
L4.1a-D  control the dynamic far tail, gap >= m(Z);
L4.1b    control the transition gaps 1 <= gap < m(Z);
L4.2     control the local gap-zero block.
```

`L4.1a-D` is now proved with explicit cutoff-independent hard-shell constants
and a moving-gap Young remainder linear in `Z`; the finite numerics remain
regressions rather than part of that proof. `L4.1b` and `L4.2` remain open.

The far-tail objective is now parameterized by its minimum dyadic gap. For
`gap >= 3` (`high/low > 8`), the absolute stopping tolerances in the optimizer
were invalid: the fourth-power objective can be nonzero far below `1e-30`.
Replacing them with relative machine-epsilon tests exposed a genuine gradient
branch. At `E(0)=1`, initial `H4^2<=100`, `nu=0.1`, and `T=0.01`, the validated
branch reaches `J_gap>=3=1.29741e-20` at K7. Its K6 projection retains `81.1%`;
zero-padding the earlier K6 winner to K7 retains `98.8%`. Zero-padding the K7
winner to K8 retains `99.9775%` in the fine-time-step check. Thus the branch persists
after the gap-three band opens, although this finite calculation neither finds
a global maximum nor proves a decay exponent.

The separately optimized `gap>=2` value is `3.21556e-11`, about `2.48e9` times
the current `gap>=3` value. Because the objective contains the fourth power of
the signed tail stretching, this ratio cannot be inserted directly as a
paraproduct decay exponent. The analytical target remains a conventional,
cutoff-independent estimate of the form

```text
sup_N integral_0^T |V_{N,gap>=m}|^4/(Z_N P_N^3) dt
    <= C(E(0), nu, T, u_0) 2^(-alpha m),  alpha > 0,
```

derived without assuming a uniform high-Sobolev solution bound. The fixed-gap
envelope alone does not close in time, but it is now one input to the proved
dynamic-tail split below.

One exact factor in that envelope is now certified. For `k=p+q` with `p` the
unique low advecting wave, define

```text
T(p,q,k) = Re <u_k, i(q dot u_p)u_q>.
```

Reality and incompressibility pair `(p,q,k)` with `(p,-k,-q)` so that

```text
T(p,q,k)+T(p,-k,-q)=0,
|k|^2 T(p,q,k)+|q|^2 T(p,-k,-q)
    = (|k|^2-|q|^2)T(p,q,k),
abs(|k|^2-|q|^2) <= |p|(|k|+|q|).
```

Thus the low-advecting gap-`m` block has an exact one-factor `2^-m`
commutator gain relative to its raw two-derivative weight. The low-advected
role has its derivative on the unique low wave, while the low-target role has
the low enstrophy weight and gains two low/high factors. `TriadTailEnvelope`
reconstructs the full separated signed ledger from these three roles with
`3.21e-19` relative error in self-test and verifies every termwise amplitude
and normalized frequency inequality with ratio at most one.

The resulting hard-shell estimate is

```text
|V_gap>=m| <= (C1 2^(-m/2)+C3 2^(-3m/2)) Z^(3/4) P^(3/4),

C1=192(2+2^(5/2))sqrt(2),
C3=64(2+2^(5/2))sqrt(8/7).                                  (FT-1)
```

Its fixed-gap fourth-power density is bounded by `C 2^(-2m) Z^2`, which is
not time-integrable from the energy inequality. Candidate F004 records this
obstruction. Applying Young before taking the fourth power instead gives

```text
|V_gap>=m| <= (nu/4) P + C nu^(-3) 2^(-2m) Z^3.
```

The state-dependent split

```text
m(t) = m0 + ceil(log2(max(1,Z(t))))
```

reduces the cubic remainder to `C nu^(-3) 2^(-2m0) Z(t)` without assuming a
future bound for `Z`. `MovingGapController` checks the integer inequality and
the scaling certificate checks the exact Young exponents. `DyadicShellBounds`
separately certifies the scalar sequence sums with constants `sqrt(2)`,
`sqrt(8/7)`, and one for the interpolated high moment.
`PeriodicShellGeometry` supplies the Fourier lattice count, adjacent-shell,
and role constants; `FarTailClosure` verifies the explicit Young remainder.
This closes the dynamically selected far tail, uniformly in the Galerkin
cutoff. The active L4 task is now the remaining local/transition block of
`O(log Z(t))` dyadic gaps. Plain band counting is rejected as F005 because it
leaves `Z^3 log(1+Z)^4`. See
`proof/l4/lemmas/dynamic-far-tail/README.md`.
