# Proof roadmap and restart point

Last updated: 2026-08-02

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

The current restart point is the **frozen-data local L4 trajectory lemma**.
For one initial datum set `k0=sqrt(Z(0)/E(0))`, `B0=E(0)P(0)`, and prove a
cutoff- and time-uniform upper bound for

```text
R_local(t)=4S(t)^3 Z(t)P(t)(K(t)+G(t))
 / [k0(S(t)^4 Z(t)^2P(t)+B0 Z(t)^3P(t)^4)].
```

The exact maximum-on-trajectory adjoint is implemented and validated. Warm
continuation from K2 gives the current K3 lower bound `8.53498799310e-4` on
`[0,0.5]`, attained at `t=0.298`; dt halving changes it by `8.16e-16`
relatively. This corrects an earlier secondary branch but is still a pattern
and falsification result, not the required analytic bound. Resume from the
cyclic response hierarchy and the `(m,m,2m)`/remainder/mixed decomposition in
`proof/l4/lemmas/shifted-local-density/SIGNATURE_BLOCK.md`.

Six accepted K4 continuation steps improve that value to
`8.53527437357e-4`, again at `t=0.298`, with dt-halving error `9.00e-16`.
The new hard-shell-four energy is only `1.97e-9`; the K4-to-K3 state residual
is `2.96e-4`, and the final projected objective gradient is `5.17e-5`.
This is the current finite lower bound, not a global optimum or a lemma.

The current K3 winner is not an unstructured high-dimensional state. Sixteen
quadratic response orders plus the explicit transverse `(2,1,1)` and two
oriented `(3,1,0)` orbit directions capture `99.9998628%` of its energy. The
projected state retains `99.96181%` of the trajectory objective and the same
peak time. The next proof step is to replace this finite projection fact by a
cutoff-uniform summability estimate for the response coefficients and orbit
remainders; finite-dimensional agreement alone is not the lemma.

At K4, the corresponding 19 directions capture `99.9997617%` of state energy
and `99.9331554%` of the full trajectory objective. Direct-triad and FFT
evaluation agree exactly at the serialized precision. A new cutoff-diagonal
ledger avoids comparing response orders after they hit the Galerkin wall: at
cutoff K it retains only orders `0,...,K`. For the same K4 state zero-padded
through K8, the weighted coefficient sum
`A_1.25(K)=sum_{n<=K}1.25^n |<u,b_n>|/sqrt(E)` remains between `1.12012` and
`1.15916` from K2 through K8. This is evidence for the explicit weighted
response target in
`proof/l4/lemmas/shifted-local-density/RESPONSE_DIAGONAL.md`.

The exact direct-triad tensor has now separated the response-space derivative
loss from cutoff growth. Same-radius weights at `r=1.25` grow from `0.9723` at
K2 to `5.2028` at K8 and cannot close the algebra. With input radius `R=2` and
output radius `r=23/20`, however, every ordered response interaction through
K12 obeys the same sharp finite bound

```text
sum_m (23/20)^m |<b_m,B(b_i,b_j)>|
    <= [23/(20 sqrt(3))] 2^(i+j),
```

with equality at the axis pair. The lowest open analytical subtask is now to
prove this explicit scalar inequality for arbitrary response order.

The transverse space is now graded before orthogonalization and extended by
the largest missing bilinear products. With all three explicit orbit families
and sixteen closure directions, the invariant degree-block projected
constants at K3/K4/K5 are `1.04467`, `1.04894`, and `1.04894`; adding the
separately maximized shell complement gives finite bounds `1.17088`, `1.30595`,
and `1.32010`. The lowest open complement subtask is to prove a uniform
dimension/support count for this indefinitely generated graded orbit tree and
sum it in the `R=2 -> r=23/20` gap. A strict response-degree support shortcut
has been rejected: K5 contains a coefficient of magnitude `0.409773` two
degrees above `d_left+d_right+1`. The proof therefore also needs quantitative
off-diagonal decay of the response tensor, not just a count of nominal support.
See
`proof/l4/lemmas/shifted-local-density/RESPONSE_TENSOR.md`.

The response route also has a newly certified limitation. At radius three,
the measured diagonal grows from `A_3(6)=2.14594` to `A_3(8)=5.97485`, and
the eighth response vector keeps `95.01%` of its energy in the lower half of
the physical shells. Response order is therefore not a proxy for Fourier
frequency, so smooth initial data alone does not provide the exponential
response weight. This route remains a structural extremizer model rather than
the current universal closure.

The universal dominant-block calculation has progressed independently. The
complete `(m,m,2m)` family is an equal-length orthogonal incidence graph.
Target-wise degree `O(R)` gives `||B_d||_2 <= C R^(3/2)E`, and every one-shell
entry of the closed `K_d+G_d` bracket has scale `R^5E^2`, gaining one half
derivative over LQC-3. The naive attempt to absorb the global normalization
sum into the structural shell sum is false: two shells with high-shell energy
`L^(-11/4)` make that intermediate ratio grow as `L^(1/8)`. Direct estimates
avoid the obstruction: `S <= C Z^(5/4)P^(1/4)` and
`T <= C Z^(1/4)P^(5/4)`, so both normalization entries are bounded by
`C Z^(3/2)P^(1/2) <= C Z^(5/4)P^(3/4)`. Together with neighbor-shell locality
of the structural terms, this proves the complete doubling-family block
cutoff-independently. Fixed-angle incidence also proves the complete
`(m,m,3m)` family, and the general plane--sphere argument proves every fixed
primitive projective ray `(am,bm,cm)`. The unresolved issue is a uniform or
summable estimate over all primitive shapes. The remainder has also been
reduced algebraically. Writing `W=B_rem-cAu` and
`D=B_rem(u,Au)-[x->B_rem(x,u)]^*Au`, two exact square completions give

```text
K_rem+G_rem
=-||A^(1/2)(W-(1/2)A^(-1)D)||^2
 +S^2/(2Z)+c^2H3+c<Au,D>+(1/4)||A^(-1/2)D||^2.
```

The direct VJP also certifies `D=-[dB_rem(u,u)]^*Au`. A growing dense
absolute-value branch is negative because of the displayed square. The first
signed-LQC3 continuation appeared to stabilize at `0.102959`, but a
cross-objective warm start finds a dense positive branch growing from
`0.103623` at K4 to `0.175398` at K6. Its normalized stretching is only
`O(1e-6)`, making the actual local source `O(1e-20)`. Maximizing only the
positive envelope and retaining a fixed fraction of the negative square also
produce growing branches. Therefore the lowest universal analytical subtask
is the exact joint bracket--shape tradeoff, not an independent bound on D or
the signed LQC3 quotient. Since
`R_rem=[(K_rem+G_rem)S_full/(Z^2P^2)] 4x^2/(1+x^4)` and the last factor is at
most `2`, the current remainder lemma is the power-one estimate
`|(K_rem+G_rem)S_full| <= C Z^2P^2`. The mixed block follows. See
`proof/l4/lemmas/shifted-local-density/DOUBLING_QUARTET.md` and
`proof/l4/lemmas/shifted-local-density/PROJECTIVE_QUARTET.md` and
`proof/l4/lemmas/shifted-local-density/REMAINDER_QUARTET.md`.

The first explicit L4 candidate in the reduction was

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
`<gradient C_local,RHS>` in about `0.22 s` at K6; the later complete role
ledger takes `0.48 s` on twelve threads. Direct and FFT RHS backends agree in
the main serialized values and within `1.1e-19` in the expanded budgets. Its
normalized rates on the four optimized
K3--K6 winners are `7.2396e-4, 7.5341e-4, 7.7546e-4, 7.8151e-4`; the K6 value
is only `6.10e-4` above the short-horizon average relatively. This removes the
need to use RK4 for first-pass SLD falsification. The unresolved work is still
analytical: split this exact derivative into local signature/triad terms and
prove a cutoff-independent upper bound, or produce a scalable counterfamily.

That first derivative split is now exact and implemented. With signed local
stretching `S`, the engine reconstructs

```text
C'_local = 4 S^3 S'/(Z P^3)
           - S^4 Z'/(Z^2 P^3)
           - 3 S^4 P'/(Z P^4)
```

from independently differentiated `S`, `Z`, and `P`; its K6 relative residual
against `gradient C dot RHS` is `1.18e-18`. It also certifies the exact PDE
identities `Z'_NL=-2S_global`, `Z'_nu=-2nu P`, and `P'_nu=-2nu H3`. The
stretching derivative is split into the outer, advecting, and advected
Frechet slots. On the K6 extremizer they are `-0.609617`, `+0.016053`, and
`+0.328901` for the nonlinear direction, so a large cancellation would be
lost by separate absolute-value estimates.
The pure local part of the outer slot is exactly
`-||A^(1/2)B_local||_2^2=-0.603413126167`; the new quartic identity ledger
certifies this with zero serialized residual rather than treating its sign as
numerical evidence.

After normalization by `(C_local+B0) k0 Z`, the K6 budget is `6.87169e-4`
from the local nonlinear RHS, `6.48121e-5` from the nonlocal nonlinear RHS,
and `2.95310e-5` from viscosity. The local block therefore supplies `91.38%`
of the nonlinear source on the current extremizer. A same-state K6--K8 audit
changes only the nonlinear stretching channel and converges from
`7.8151205976e-4` to `7.8151208069e-4`.

The active restart point is no longer differentiation. It is the complete-
quartet symmetrization of the local-local contributions to `S'` together with
the nonlinear `P'` term in the denominator-free form

```text
4 S^3 S' Z P - S^4 Z' P - 3 S^4 Z P'
 <= A k0 (S^4 Z^2 P + B0 Z^3 P^4).                 (SLD-1P)
```

The proof must preserve the cross-slot cancellation recorded by the ledger;
a separate norm bound for every slot is not yet a closure because it discards
the signed quartet structure.
On the K6 oracle, the outer negative-square source is `2.30680e-3` after SLD
normalization, while the advecting, advected, enstrophy, and palinstrophy terms
reduce it to `6.87169e-4`, a `70.21%` cancellation. The denominator-free
SLD-1P evaluation agrees with the shifted-log route to `4.74e-19` relatively,
so future quartet algebra can be regression-tested without differentiating a
quotient.

The local quartet is now resolved by target mode, exact eigenshell, hard
shell, and dyadic annulus. On the optimized K6 state, `70.66%` of component
cancellation occurs within one dyadic annulus and essentially none occurs
between annuli. Control states show the opposite placement, so this is routing
information rather than a universal cancellation lemma.

One infinite-dimensional estimate has nevertheless been completed. If
`B_L=B_local(u,u)`, elementary local-shell geometry and lattice counting give

```text
||(A^(1/2)B_L)_j||_2^2 <= 46656 R_j^7 (E_j^near)^2.
```

Three-shell overlap and interpolation
`sum R_j^3 E_j <= sqrt(ZP)` then give the cutoff-independent global lemma

```text
||A^(1/2)B_L||_2^2 <= 7264120.5 sqrt(Z) P^(3/2).       (LQE-2)
```

This removes the former high-Sobolev concern for the isolated outer slot.
The code certifies every geometry and summation stage; the conventional proof
is in `proof/l4/lemmas/shifted-local-density/SHELL_ENVELOPE.md`.

The outer and advected slots are also no longer independent. Transport
skew-symmetry gives

```text
-<A B_L,B_L>-<A u,B_L(u,B_L)>
 =-<B_L,A B_L-B_L(u,A u)>,
```

whose Fourier multiplier satisfies
`abs(|k|^2-|q|^2)<=|p|(|k|+|q|)`. The K6 identity residual is below `2e-19`.
The active restart point is therefore narrower: bound the remaining
advecting-slot and nonlinear-palinstrophy combination together with this
commutator in SLD-1P. The completed outer estimate must be reused rather than
rederived with a higher solution norm.

An alternative exact grouping has now absorbed both denominator derivatives.
For `h=-B_L`, their combination with the outer slot is

```text
<A h,B_L-Su/(2Z)-3S A u/(2P)>
 =-||A^(1/2)(B_L-3S A u/(4P))||_2^2
   +S^2/(2Z)+9S^2H3/(16P^2).
```

After the common scalar coefficient is restored, this projected pairing plus
the advecting and advected slots reconstructs the entire local quartet. The
K6 normalized residual is `3.85e-18`. This is the current minimal algebraic
target. Because the common coefficient contains signed `S^3`, neither the
negative square nor its positive remainders may be discarded uniformly. The
next lemma must be a sign-aware estimate for this three-entry combination,
not another five-way triangle inequality.

Merging that grouping with the outer--advected commutator leaves only two
entries. Define

```text
K=-<B_L,A B_L-B_L(u,A u)>+S^2/(2Z)+3S<A B_L,A u>/(2P),
G=<A u,B_L(-B_L,u)>.
```

The exact local polynomial numerator is `4S^3ZP(K+G)`. At K6 it equals
`0.00264627230688`; the independent five-entry calculation agrees to
`3.20e-19` relatively. The active local lemma is now SLD-1P-L in
`proof/l4/lemmas/shifted-local-density/REDUCED_QUARTET.md`: prove a
cutoff-independent upper bound for this signed two-entry expression. This is
the next restart point after any failed computational or analytical attempt.

The exact Young reduction first identified the following sufficient absolute
estimate:

```text
|K+G| <= C k0 B0^(1/4) Z^(5/4)P^(3/4).               (LQC-3)
```

Multiplication by `4|S|^3ZP` and weighted AM--GM imply SLD-1P-L with
`A_local=3C`. Compared with the elementary quartic scale
`Z^(1/2)P^(3/2)`, LQC-3 requires a `(Z/P)^(3/4)` depletion supplied by the
signed two-entry structure. The K6 finite-state ratio is `0.127982322632` and
the algebraic identity error is `8.00e-20`; neither number proves uniformity.
`proof/l4/lemmas/shifted-local-density/CLOSURE_TARGET.md` is the exact current
algebraic target, but the absolute estimate is stronger than the signed local
polynomial inequality.

For the closed doubling remainder, an exact-gradient absolute LQC-3 search
finds a dense negative branch: its magnitude grows from `0.05903` at K2 to
`0.42906` at K5. A signature-row ledger reconstructs the K5 bracket at
`8e-19` relative error and measures 549.49 effective contributors with
`0.999840` same-sign alignment. This is not a counterexample to the required
upper bound. The first signed LQC-3 search padded a `0.102959` K5 low-mode
winner through K8, but a later cross-objective warm start rejects that basin
as a worst-case proxy: it gives `0.103623`, `0.136748`, and `0.175398` at
K4--K6 with 452.20 effective K5 signature contributors. This stronger branch
has nearly zero stretching, so its actual local SLD ratio is below `6e-20`.
Direct optimization of the exact remainder block product instead stays near
`0.00022068` from K3 through K6. The bracket--shape factorization in
`proof/l4/lemmas/shifted-local-density/REMAINDER_QUARTET.md` is therefore the
current remainder restart point. Its exact scalar reduction leaves the
energy-independent power-one tradeoff
`|(K_rem+G_rem)S_full| <= C Z^2P^2`; the corresponding exact-gradient roots
remain between `0.00601` and `0.006365` through K7. Projective grouping shows
that `95.42%` of the K7 absolute total lies on one primitive ray. After the
proved doubling and triple families are removed, the K5 tail has only `3.027`
effective projective shapes. Every fixed ray and the projective square
function are proved cutoff-independently. A coherent fan proves that a
standalone synthesis estimate cannot be uniform: its ratio is at least
`K^2/640`, while its stretching is exactly zero. The remaining analytic gap
is therefore the stretching-aware cross-ray tradeoff in RQ-11, not unsigned
uniform summation.
The exact self/cross quartet ledger shows that the K7 power-one winner splits
as `0.00614715` same-ray and `0.00021762` unequal-ray. A new cross-only exact
gradient objective stays near `0.0011166` from K4 through a K8 lift, while its
K5 winner exhibits opposite-sign self/cross cancellation. A per-ray cross
attribution has only `2.103` effective rays on that K7 state and `3.583` on
the K6 cross-only winner. The plane--sphere incidence proof now applies to
the union of any fixed finite ray set, so its complete internal self+cross
quartet closes cutoff-independently with gain `-1/2`; the replayable 13-ray
certificate is documented in `FINITE_PROJECTIVE_FAMILY.md`. The active
analytic subtask is no longer the internal cross term of a fixed core. It is
the core--growing-tail and tail--tail part of PCQ-3, with a constant uniform
under movement of the core boundary.

That boundary is now a first-class exact-gradient objective rather than a
residual inferred from unrelated winners. For the canonical finite core
`F_H={primitive feasible (a,b,c): max(a,b,c)<=H}`, the engine subtracts its
complete internal quartet and maximizes only
`|(J_FT+J_T)S_full|^2/(Z^4P^4)`. The gradient error is `8.94e-12`. The H=8
branch stabilizes near `0.0010076` by K8. Separately adapted K8 branches at
H=16,32,64 reach `0.0008591`, `0.0006330`, and `0.0003079`, giving a finite
lower-branch height slope `-0.557`. Because the optimizer relocates toward
higher primitive shapes, this decay cannot yet be used as an upper bound.
The growing boundary now has an exact dyadic primitive-height matrix. Its
absolute five-component envelope obeys the proved finite Schur reduction
`sum_(i<=j)e_ij <= R sum_i e_i`. On the principal K8 stress states the
normalized row sum is `2.89`--`3.85`. Separately optimized diagonal states
give signed and absolute-envelope finite height slopes `-0.374` and `-0.435`.
The lowest open analytical subtask is now precise: prove a cutoff-uniform
bound for `R` and a cutoff-uniform summation of the diagonal envelopes. The
finite slopes are only lower-branch evidence, not either upper bound. See
`proof/l4/lemmas/shifted-local-density/DYADIC_PROJECTIVE_HEIGHT_SCHUR.md`.

The absolute closure ratio was subsequently reverse-differentiated and
optimized directly. A K1 branch reaches `C_state=1/3`, but its stretching is
only `1.23e-11`; the `S^3` factor makes it irrelevant to the actual local
polynomial source. The engine therefore now optimizes the signed quotient

```text
R_local=4S^3ZP(K+G)/[k0(S^4Z^2P+B0Z^3P^4)]          (LQC-5)
```

instead of treating LQC-3 as the only search objective. Its analytic gradient
agrees with central differences to `6.77e-12`. Twelve-start L-BFGS continuation
gives `R_local=7.60409e-4` at K1 and `7.95364960e-4` at K8; gains after K3 are
below `4e-10`, and the final adjacent projection residual is `1.99e-6`.

The K8 extremizer is a low-shell, predominantly heterochiral state. Six axis
modes and twelve face-diagonal modes dominate, while the `(1,1,2)` local
signature supplies `98.8%` of its signed transfer. The explicit two-basis
cyclic-shear/quadratic-response ansatz reaches `7.75010870e-4`, or `97.44%` of
the full value. The active restart point is to prove a signed bound for this
dominant block plus a cutoff-uniform estimate of the orthogonal shell
remainder. The present direct quotient uses the evaluated state as the datum,
so a trajectory proof must additionally freeze `k0,B0` at time zero. These
finite optimizations do not complete L4 or the Clay problem.

The dominant block is now an exact operator decomposition rather than a table
classification. `L_d` selects every local squared-frequency signature
`(m,m,2m)` and `L_r=L_local-L_d`. Independent evaluation of all mixed
bilinear terms certifies

```text
K+G = F_d + F_r + F_cross
```

to `2.20e-18` relatively on the K8 winner. Their direct SLD contributions are
`7.70182e-4`, `1.21606e-5`, and `1.30224e-5`. The next proof attempt must
derive a scale-uniform signed estimate for `F_d` and a summable joint estimate
for `F_r+F_cross`; separate triangle inequalities are not accepted because
they restore the known frequency loss.

The tempting global sign shortcut has now been eliminated. Exact-gradient
search gives a positive normalized bracket `c=0.0740187069851` for the full
local operator, and positive examples also exist for both sides of the
doubling-family split. Those states have nearly zero stretching, so the
correct restart point is the joint factorization `R_j=c_j phi(x)`, not the
sign of `c_j` alone. The common-normalization block gradients pass central
differences at `1.29e-11` or better. Separate cutoff scans stabilize near
`7.72930e-4` (doubling family, K6), `2.20683e-4` (remainder, K6), and
`1.34936e-4` (mixed, K6). The remainder and mixed blocks therefore need a
real joint estimate even though they are small on the original K8 branch.

The static absolute search has also isolated a sharp candidate rather than
only a decimal plateau. The normalized six-mode cyclic axis shear satisfies
`E=Z=P=1`, `S=0`, and `K+G=-1/3` with `1.08e-19` identity error. Its padded
branch stays at `|c|=1/3` from K1 through K8 under 12-start searches. The
candidate inequality

```text
|K+G| E^(1/4) <= (1/3) Z^(7/4)P
```

is therefore sharp if true, but remains unproved.

More importantly, the engine now removes the `t=0` shortcut. A separate
trajectory adjoint freezes `k0=sqrt(Z(0)/E(0))` and `B0=E(0)P(0)`, evolves
with RK4, and includes the direct initial-data derivatives of both constants
in the reverse pass. Its gradient error is `5.62e-12`, and the zero-step
limit matches the static gradient to `1.52e-19`. At K3 and `nu=0.1`, optimized
terminal ratios rise from `7.98918e-4` at `T=0.01` to `8.46863e-4` at
`T=0.20`, with dt-halving errors below `4.3e-16`. Thus the current restart
point is a uniform-in-time bound for the frozen-data joint quotient, not the
static `1/3` conjecture alone.

The maximum-on-trajectory adjoint removes endpoint bias. On the current K3
search over `0<=t<=0.5`, it selects `t=0.2175` and returns the refined lower
bound `8.48675785e-4`; coarse/refined checkpoint maxima differ by `3.10e-6`
relatively. The optimized initial state remains low-shell: `99.9723%` of its
energy lies in the first hard shell and the `(1,1,2)` signature supplies
`99.42%` of the coherent signed transfer. The nonzero projected gradient
`3.55e-4` prevents treating this number as a converged maximum.
Evolving that initial state to the refined peak and retaining its initial
`k0,B0` splits the quotient into `8.18560410e-4` (doubling family),
`1.52599847e-5` (closed remainder), and `1.48553902e-5` (mixed), with
`9.98e-19` reconstruction error. The family therefore remains the primary
trajectory lemma, while the latter two terms require one joint summable
estimate.
