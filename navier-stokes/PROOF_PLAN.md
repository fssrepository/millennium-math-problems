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

The local block is now decomposed into all eight helical sign sectors. On the
K3 local-objective endpoint, `V_local=-0.227783`, while its homochiral part is
only `-3.71e-6`; the heterochiral sectors carry the signed extremum. This is
not an exact identity: pure positive- and negative-helicity states have
nonzero homochiral local stretching, so F006 rejects universal homochiral
cancellation. The next computational target is a sector-selective local
objective; its exact static and checkpointed trajectory gradients now pass at
`6.71e-12` and `3.56e-13`. Twelve parallel restarts reach static critical
density `1.33682e-4` and short-time critical integral `5.07257e-8`.

The same-state cutoff continuation is now complete through K6. At
`nu=0.1`, `dt=0.001`, and `T=0.01`, its heterochiral local critical integral is

```text
K:       3                  4                  5                  6
J: 4.39045255695e-6   4.39045682754e-6   4.39045682741e-6   4.39045682741e-6
```

The K5-to-K6 relative difference is `2.67e-15`. Twelve large-perturbation
static multistarts through K5 also select the smooth original continuation,
with critical density saturating near `4.40316e-4`. Thus the current branch is
cutoff-stable and does not falsify the dynamic lemma. This finite convergence
does not provide the uniform analytical estimate.

The next analytical target is a quantitative heterochiral depletion statement
strong enough to supply the missing `Z^(-1/2)` factor. The next computational
task is no longer another lift of this same branch: it is to encode and
falsify candidate dynamic mechanisms (sector cancellation, relative-helicity
coercivity, or time-integrated shell transfer) that could produce that factor.
See `proof/l4/lemmas/local-helicity/README.md`.

The first two mechanisms are now excluded on the converged branch: it has
balanced positive/negative helical energy and a nonvanishing heterochiral
signed-to-absolute ratio near `0.202`. `HelicalGapLedger` places
`-0.226934` of signed heterochiral transfer in gap zero, `3.45e-4` in gap one,
and only `-4.61e-9` in gap two.

For the remaining gap-zero block, complete-triad symmetrization proves the
cutoff-independent identity

```text
|V_triad| <= (max |k_i|^2-min |k_i|^2) sum_i |T_i|.          (LS-1)
```

Equal-frequency triads therefore contribute exactly zero. The K3 endpoint,
however, places `-0.207019` of its total `-0.226938` local transfer in the
squared-frequency-spread interval `(1/4,1/2]`. Consequently LS-1 is a genuine
partial lemma but supplies no small scale factor on the dominant broad-spread
class. The restart point is now a dynamic/time-integrated estimate for those
broad-spread gap-zero triads. See
`proof/l4/lemmas/local-triad-symmetry/README.md`.

The objective/adjoint now masks this exact broad class. A 12-restart K3 search
reaches `4.39098e-6`; the narrow and equal-frequency objectives on the same
state are only `2.71e-27` and `8.01e-89`. The broad same-state K3--K6 scan has
last relative difference `2.66e-15`. This removes cutoff and optimizer
artifacts from the current branch and isolates the analytical gap still
further: prove a trajectory-integrated estimate specifically for broad-spread,
gap-zero heterochiral triads.

Signature aggregation makes that target concrete. The squared-length family
`(1,1,2)` carries `91.27%` of the broad K3 endpoint's signed local transfer,
and `(1,1,2)` plus `(2,2,6)` carry `99.74%`. The first family is the geometry
of equal-length orthogonal input waves. The next lemma attempt must therefore
bound the time-integrated coherent multiplicity of scaled orthogonal local
triads; an estimate for one isolated triad would miss the lattice-count growth.

That exact signature family is now controlled analytically. For fixed `p`, the
integer constraints `|p|=|q|` and `p dot q=0` admit at most `2(2K+1)` partners
in a cutoff-K cube; the same bound holds for a fixed target. The resulting
hypergraph Cauchy bound is

```text
|V_(r,r,2r),R| <= C R^(7/2) E_R^(3/2)
                 <= C R^(-1/2) E(0)^(1/2) P_R.
```

Hence viscosity absorbs this family on all sufficiently high shells, and only
finitely many low shells remain. This is a proved partial local lemma, not a
finite numerical inference. It covers the signature carrying `91.27%` of the
current K3 endpoint, but that experimental fraction cannot remove the other
signatures. The restart point is the degree/multiplicity classification of
nearby broad-spread squared-length triples. See
`proof/l4/lemmas/local-orthogonal-triads/README.md`.

The exact multiplicity threshold is now encoded. If the coherent partner
degree is `D(R)=O(R^d)`, the transfer has frequency power `3+d/2`.
Unrestricted local interactions have `d=3` and power `9/2`; viscosity has
power `4`; the critical degree is `d=2`; the controlled orthogonal family has
`d=1` and power `7/2`. A closing extension must therefore produce effective
degree `d<2` on the unresolved broad class.

The extension to all fixed squared-length signatures is now proved. Every
ordered `(r,s,t)` family has input and target degree at most `2(2K+1)`. More
importantly, summing the squared per-signature estimates before Cauchy gives

```text
(sum_sigma |V_sigma|^2)^(1/2) <= C K^(7/2) E_K^(3/2).
```

With
`N_eff=(sum_sigma |V_sigma|)^2/sum_sigma |V_sigma|^2`, the entire local block
is bounded by `C sqrt(N_eff) K^(7/2) E_K^(3/2)`. Therefore
`N_eff=O(K^mu)` closes exactly when `mu<1`. Replay of the broad K3 endpoint
gives `N_eff=1.1902068038` and dominant coherent fraction `0.912689851446`,
despite 106 signatures above the numerical threshold. This is a sharper
restart point than raw signature counting, but the unrestricted pointwise
claim has now been rejected. A 12-worker flat-spectrum search reaches
`N_eff=1748.62` at K6 and has finite-range slope `4.223`; energy and
incompressibility alone do not force sublinear participation. The viable
target retains inter-signature signs. With

```text
A_sig=|sum V_sigma|/(sum |V_sigma|^2)^(1/2),
```

the square-summed lemma closes whenever `A_sig=O(K^alpha)`, `alpha<1/2`.
The same K1--K6 search gives K6 maxima `2.33549` (flat) and `2.73471`
(outer-half-flat), but this random screen is not adversarial enough. The exact
analytic-gradient search raises the flat values to `4.59, 10.09, 19.65,
32.01, 49.39` on K2--K6, with exponent `2.168`, and rejects pointwise LSF-4.

The remaining identity preserves the coupled magnitude:

```text
V^4/(Z P^3) = A_sig^4 (sum_sigma |V_sigma|^2)^2/(Z P^3).
```

On the broad same-state K3--K6 trajectory, its factorization residual is below
`3.1e-19`, its critical integral is `4.39122e-6`, and the last cutoff
difference is `1.51e-15`; halving the RK4 step changes it by only `4.44e-7`
relatively. The analytical restart point is a
trajectory-integrated bound for this product tied to a fixed smooth initial
datum; neither factor has a sufficient pointwise bound by itself. See
`proof/l4/lemmas/local-signature-families/README.md`.

The coupled local objective has now also passed a true independent dynamic
multistart screen. Twelve concurrent exact-adjoint starts with eight L-BFGS
steps at each cutoff give

```text
K:          3             4             5             6
J_local: 4.392286e-6  4.392821e-6  4.392876e-6  4.392883e-6
top E:   2.475e-7     2.130e-8     3.557e-11    1.133e-13
```

The smooth continuation wins every cutoff; at K6 an independent competitor
is only `2.91e-5` lower relatively. This removes the former single-dynamic-
start weakness but remains finite evidence. The restart point stays the same:
prove or falsify a cutoff-uniform trajectory estimate for the coupled product,
next under longer time horizons and lower viscosity. See
`proof/l4/analysis/local-signature-coupled-integral/README.md`.

The first longer-horizon local multistart at K6 gives
`J_local(0.02)/(2J_local(0.01))=0.9974002` and normalized log-Q gain
`0.9655509`. Thus the current branch accumulates the critical integral
slightly slower than linearly; no accelerating short-time branch appears.
Lower-viscosity continuation remains the next finite falsification check.

At the same K6 horizon, reducing viscosity from `0.1` to `0.02` changes the
optimized local integral by a factor `1.0038754`, while the log-Q gain grows by
`1.3093501`. Twelve independent starts again select the smooth continuation,
with the next branch only `2.43e-5` lower. This still does not prove the
uniform lemma, but neither time nor viscosity continuation has produced a
short-time concentration obstruction. The next diagnostic is the timewise
correlation inside the exact coupled signature factorization.

That correlation is nearly perfect on the optimized branch:
`corr(log A_sig^4, log(R^2/(ZP^3)))` lies between `-0.999980` and
`-0.999978` on its K3--K6 trajectories, and is `-0.999848` on the
lower-viscosity K6 trajectory. However, a 240-state parallel one-step
adversary finds 41 states in which both factors grow simultaneously. The same
counts survive halving the RK4 probe step. F009 therefore rejects universal
factor antimonotonicity. Any viable proof must derive a trajectory-integrated
estimate tied to the fixed smooth datum, not a pointwise sign identity.

The coupled density is not pointwise monotone either. An exact discrete
endpoint adjoint for

```text
C_local(u(T))-C_local(u(0)),
C_local=|V_local|^4/(ZP^3),
```

passes its centered-difference test at `5.23e-12`. Twelve-start projected
L-BFGS searches at `T=0.001` find positive increases
`9.64805e-7, 1.09941e-6, 1.12519e-6, 1.13495e-6` on K3--K6. Time-step errors
are below `2.5e-14`; top-shell energy decays with fitted exponent `-8.18`, and
the K5-to-K6 projection residual is `1.74e-3`. F010 therefore rejects the
nonincreasing-density route. The increase is already flattening in cutoff, so
the data identify a real transient-growth mechanism but not a singularity or
a failure of the time-integrated lemma. The next machine task is horizon
continuation of the refined K6 branch; the analytical restart point remains a
cutoff-uniform integral estimate rather than a pointwise sign argument.

The first horizon continuation resolves that diagnostic. At K6,
`Delta C_local(0.002)/(2 Delta C_local(0.001))=0.999368654`, so the positive
short-time growth is nearly linear rather than accelerating. The signature
factor correlation on this endpoint-optimized trajectory is `+0.999994857`,
with both factors peaking at the final sample and factorization residual
`2.32e-19`. This shows that the strong anticorrelation selected by the
integral maximizer is branch-dependent even at substantial critical density.
The next falsification axis is lower viscosity; absent accelerated or cutoff
growth, the proof restart remains an integrated estimate exploiting more than
pointwise factor signs.

Reducing viscosity fivefold to `nu=0.02` at K6 and `T=0.002` lowers the
optimized absolute increase by a factor `0.909393069` and its logarithmic gain
by a factor `0.823480963`. The endpoint error under time-step halving is
`2.00e-15`. Thus neither horizon nor viscosity continuation reveals an
accelerating branch. The engine now differentiates the more lemma-directed
objective `log(C_local(T)/C_local(0))`; its central-difference error is
`2.08e-11`, and cutoff searches report the normalized average
`log(C(T)/C(0))/(T k0 Z(0))` directly.

The unshifted objective immediately exposes a zero-set obstruction rather than
a high-density extremizer. It drives `C_local(0)` to the `1e-30` numerical
guard and reaches log gain `31.81` at K6, with projected gradient `3.8e7`.
Consequently a purely multiplicative estimate is not a stable analytical
restart point near vanishing local transfer. The exact adjoint and CLI now
support an additive density shift `B` and optimize
`log((C_local(T)+B)/(C_local(0)+B))`; the shifted gradient passes centered
differences at `2.79e-10`. With `B=1e-4`, twelve-start K3--K6 searches give
normalized shifted gains `3.889, 4.068, 4.119, 4.135`. Top-shell energy has
fitted exponent `-7.48`, and the K5-to-K6 projection residual is `4.14e-3`.
This is a smooth finite candidate constant, not a uniform estimate. The active
diagnostic is its dependence on `B`; only a shift tied to already controlled
initial data and correct Navier--Stokes scaling could enter a proof.

The first such explicit candidate is now encoded. Set
`B0=E(0)P(0)` and `k0=sqrt(Z(0)/E(0))`. Exact rational checks show that
`C_local` and `B0` both have amplitude degree four and Navier--Stokes scaling
exponent two, while both sides of

```text
d/dt log(C_local+B0) <= A(u0,nu,T) k0 Z(t)                 (SLD-1)
```

have scaling exponent two. If SLD-1 holds uniformly in cutoff, the energy
identity gives a uniform pointwise bound for `C_local`, hence its finite-time
integral. SLD-1 is not proved; the exact scaling and conditional Gronwall
closure only establish that it is a logically viable target. The restart
point is now an analytic differentiation of `C_local` and a triadwise bound
for its derivative with no future high-Sobolev input. See
`proof/l4/lemmas/shifted-local-density/README.md`.

The state-dependent `E(0)P(0)` adjoint passes centered differences at
`2.06e-10`. Twelve-start K3--K6 searches give normalized SLD rates
`7.240e-4, 7.533e-4, 7.751e-4, 7.810e-4`. A same-state K6--K8 zero-padding
audit changes the objective by `-7.43e-6` and then `-2.14e-9` relatively.
Thus the first exact finite adversary does not falsify SLD-1 and has converged
in cutoff on its winning branch. Further lifts of this branch are no longer
the active task; the restart point is the symbolic/local-triad expansion of
`dC_local/dt` needed to prove or reject SLD-1 analytically.

The new `shifted-density` oracle computes that instantaneous derivative as
`<gradient C_local,RHS>` in about `0.22 s` at K6. Direct and FFT RHS backends
agree in the serialized result. Its normalized rates on the four optimized
K3--K6 winners are `7.2396e-4, 7.5341e-4, 7.7546e-4, 7.8151e-4`; the K6 value
is only `6.10e-4` above the short-horizon average relatively. This removes the
need to use RK4 for first-pass SLD falsification. The unresolved work is still
analytical: split this exact derivative into local signature/triad terms and
prove a cutoff-independent upper bound, or produce a scalable counterfamily.
