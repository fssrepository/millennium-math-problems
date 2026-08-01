# Navier–Stokes Lemma Laboratory

A dependency-free C++20 research harness for the three-dimensional
incompressible Navier–Stokes equations. It combines proof-oriented algebraic
checks with a small periodic finite-difference solver on `[0, 2π]³`.

> **Scope:** finite computations cannot prove the Clay Millennium problem.
> This repository is organized around a chain of explicit lemmas: code rejects
> false candidates, verifies algebraic certificates, and records the first
> remaining analytical gap. The numerical solver is a supporting instrument,
> not the argument.

The current proof roadmap and restart point are in
[`PROOF_PLAN.md`](PROOF_PLAN.md).

## What would count as solving the problem?

The concrete target is statement **(B)** in
[Fefferman's official Clay problem description](https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf):
global smooth periodic solutions in three dimensions for every smooth,
periodic, divergence-free initial velocity and zero external force. Proving
this periodic statement alone satisfies the official problem; an additional
transfer to `R³` is not needed. The project reaches its goal only if it produces
a rigorous argument for that statement, not merely stable simulations. On the
regularity branch this requires all of the following:

1. start from arbitrary smooth, finite-energy, divergence-free initial data;
2. construct smooth Galerkin approximations and retain the standard uniform
   energy inequality;
3. prove a new cutoff-independent a priori bound in a critical regularity
   class (the open `L4` lemma in the roadmap);
4. use that bound to prevent finite-time concentration of vorticity or
   frequency flux for every finite time;
5. pass to the infinite-dimensional limit without losing the bound, then prove
   smoothness and uniqueness globally in time;
6. reconstruct a periodic smooth pressure and verify every estimate and limiting
   step under the exact assumptions of statement (B).

An alternative complete solution would rigorously construct admissible smooth
initial data whose solution develops a finite-time singularity. The repository
currently pursues the global-regularity branch. The main unresolved step is
`L4`: ordinary energy, enstrophy, and palinstrophy estimates reduce only to a
superlinear `Z³` inequality, which does not rule out blow-up. The code is meant
to isolate, falsify, and ultimately certify a genuinely stronger geometric or
dynamic flux lemma.

## Fast build

No third-party C++ libraries are required; the build needs a C++20 compiler and
CMake. The cached parallel build path is:

```bash
chmod +x build.sh
./build.sh
```

It uses all processors reported by `nproc`, retains object files for fast
incremental checks, and runs the self-tests. Override the worker count with
`NS_BUILD_JOBS=12`; set `CMAKE=/path/to/cmake` if CMake is not on `PATH`.
The equivalent manual commands are:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Engine architecture

The proof engine is split by responsibility instead of being one monolithic
source file:

- `ScalingCertificate` stores exact rational scaling results, while
  `ScalingAnalyzer` constructs the certificate;
- `SpectralGalerkin` selects direct or dealiased-FFT evaluation and controls
  kernel-level worker counts;
- `SpectralStateOps` and `SpectralStateFactory` own Fourier-state invariants,
  normalization, deterministic families, mutation, and cutoff lifting;
- `SpectralDynamics` owns direct/FFT advection, the viscous Galerkin right-hand
  side, spectral constraint projection, and RK4 time steps;
- `SpectralFftOperator` implements the dealiased forward, tangent, and adjoint
  pseudospectral kernels and is checked against direct triad summation;
- `SpectralObjective`, `SpectralAdjoint`, and `GradientAdversary` provide the
  exact `Q=D^4 Z` and `integral D^4 Z^2 dt` gradients, checkpointed reverse
  RK4 passes, and fixed-energy Riemannian gradient search;
- `InitialSobolevConstraint` projects search directions and retracts trials
  onto a cutoff-independent initial homogeneous Sobolev cap;
- `TrajectoryAnalyzer` evolves and samples trajectories, including the
  critical integral, energy balance, local/nonlocal stretching, geometry, and
  time-step refinement diagnostics;
- `TriadPartitioner` owns the single exact scale-ratio classification used by
  forward objectives, adjoints, trajectory sampling, and proof checks;
  `TriadSelection` represents arbitrary dyadic gap intervals and tails without
  adding a new hard-coded objective for every separation scale;
- `TriadLedger` resolves nonlocal stretching by dyadic gap and by whether the
  low mode is the advecting, advected, or target Fourier factor;
- `TriadCommutator` pairs low-advecting interactions with their exact
  reality/incompressibility partners and certifies the resulting low/high
  frequency gain after enstrophy weighting;
- `TriadTailEnvelope` verifies the low-advecting, low-advected, and low-target
  frequency/amplitude bounds and reconstructs the separated signed ledger;
- `HelicalTriadLedger` reconstructs local and total stretching from all eight
  curl-eigenmode sign sectors and produces pure-helicity counterexamples;
- `MovingGapController` implements the logarithmic gap choice that turns the
  post-Young far-tail remainder from `Z^3` into a linear `Z` term;
- `DyadicShellBounds` isolates the cutoff-independent scalar shell sums,
  including explicit `sqrt(2)`, `sqrt(8/7)`, and Holder constants;
- `PeriodicShellGeometry` and `PeriodicTailBound` supply explicit torus
  lattice, adjacent-shell, and final hard-shell FT-1 constants;
- `FarTailClosure` verifies the cutoff-independent moving-gap Young reduction,
  while `TransitionBlockScaling` rejects logarithmic band counting as a
  closure mechanism;
- `TriadVerifier` owns direct interaction analysis, detailed triad
  cancellation, local/nonlocal flux partitioning, and certificate aggregation;
- `StateAnalyzer` and `StateFamilyAnalyzer` measure shell decay, active modes,
  Sobolev norms, and projectivity of replayable cutoff families;
- `LemmaAdversary` schedules independent optimizer restarts;
- `ProjectiveFamily` schedules cutoff projections and switches to internal FFT
  parallelism for a single expensive high-cutoff trajectory;
- `ParallelExecutor` provides the bounded CPU worker pool shared by the two
  search drivers;
- `LemmaReporter`, `AdversaryReporter`, and `FamilyReporter` own all console
  and JSON certificate formatting through typed report objects;
- `LemmaCli` owns proof-runner argument parsing and help text.

`lemma_engine.cpp` now coordinates proof runs and keeps the integrated
self-test; state construction, trajectory diagnostics, triad verification,
forward dynamics, the discrete adjoint, constrained optimization, CLI, and
report generation are separate compilation units. The next mathematical task
is a cutoff-independent paraproduct estimate for the measured dyadic tail.
This keeps rebuilds dependency-free and makes each layer independently
replaceable.

These are ordinary non-virtual C++ classes. The numerical kernels do not use
RTTI, `std::function`, `shared_ptr`, or exceptions for control flow. Class
boundaries therefore add no intended per-mode dispatch cost; current
performance work targets state copies, temporary FFT buffers, cache locality,
and parallel work partitioning instead of removing the architecture.

## Proof-oriented lemma checks

The lemma engine checks exact scaling constraints for vortex-stretching bounds
and numerically evaluates the Fourier-triad identities used by the proof plan:

```bash
./build/navier_stokes_lab lemma \
  --triad-cutoff 3 --triad-samples 64 \
  --certificate proof/l2/l2-certificate.json
```

The first target is the family

```text
|V(u)| <= C E^a Z^b P^c,
E = ||u||₂², Z = ||∇u||₂², P = ||Δu||₂²,
V = |<Δu, (u·∇)u>|.
```

The program solves the amplitude-homogeneity and Navier–Stokes scaling
constraints with rational arithmetic. It then checks whether absorbing `P`
with Young's inequality could yield a globally closing enstrophy estimate. The
certificate records the obstruction and the Fourier cancellation residuals.

It also derives the current key-lemma target. With the scale-invariant
vortex-stretching depletion factor

```text
D = |V| / (Z^(3/4) P^(3/4)),
```

Young's inequality reduces the missing regularity step to a cutoff-uniform
bound on `integral D^4 Z^2 dt`. Reducing this further to the standard
energy-level integral of `Z` requires an exact quarter-power depletion; the
rational engine verifies that the same `1/4` requirement occurs throughout the
entire absorbable `E^a Z^b P^c` family. Establishing that dynamic depletion for
all solutions remains the open `L4` task.

## Adversarial L4 search

The static adversary mutates divergence-free Fourier states and maximizes the
strong quarter-depletion quantity `Q = D^4 Z` at fixed unit energy:

```bash
./build/navier_stokes_lab adversary \
  --cutoffs 1,2,3 --restarts 4 --generations 80 \
  --dynamic-generations 24 --dynamic-objective max-q \
  --dynamic-optimizer gradient \
  --nu 0.1 --evolve-time 0.1 --dt 0.002 \
  --certificate proof/l4/adversary/l4s-maxq-adversary.json \
  --state-dir proof/l4/states/l4s-maxq-adversary
```

Winning states are saved as Fourier coefficients, so a suspected obstruction
can be replayed or converted into an analytical family. The default dynamic
objective is the L4-A integral `integral D^4 Z^2 dt`; the exact trapezoidal
discrete adjoint projects its gradient onto the fixed-energy sphere and applies
monotone backtracking. `--gradient-method steepest|lbfgs` switches between the
normalized projected gradient and an eight-pair projected limited-memory
quasi-Newton direction. Every JSON certificate records the objective,
projected gradient norm, accepted step, Sobolev value, and line-search cost for
each iteration. `max-q` remains available for the stronger L4-S route.
`--dynamic-optimizer mutate` retains the old random search as a control, while
`hybrid` runs both. The command then directly integrates `D^4 Z^2`. Cutoff
growth attacks the strong pointwise lemma; growth of the dynamic integral
attacks the weaker time-integrated L4-A candidate.

The two terms of the local/nonlocal proof split can be attacked directly:

```bash
./build/navier_stokes_lab adversary \
  --cutoffs 3 --dynamic-generations 16 \
  --dynamic-objective critical-nonlocal-integral \
  --dynamic-optimizer gradient --gradient-method lbfgs \
  --sobolev-order 4 --sobolev-cap 100 \
  --nu 0.1 --evolve-time 0.01 --dt 0.002 \
  --threads 12 --certificate proof/l4/adversary/nonlocal.json
```

`critical-local-integral` selects the complementary objective. These masked
objectives currently use exact direct triad sums; the total objective retains
the dealiased FFT backend. Their deterministic trajectory-gradient errors are
`4.23e-12` and `3.25e-12` against central differences. The first cap sweep and
its limitations are recorded under
`proof/l4/analysis/partitioned-critical-integrals/`.

`critical-near-nonlocal-integral` selects the first separated dyadic band
`2 < high/low <= 4`; `critical-far-nonlocal-integral` selects the tail
`high/low > 4`. The latter has a nonzero K3 finite-difference gradient oracle
error of `4.84e-13`. Its first cutoff continuation is under
`proof/l4/analysis/far-nonlocal-h4-cap100/`.

Arbitrary deeper tails use the same exact direct forward/VJP kernels:

```bash
./build/navier_stokes_lab adversary --cutoffs 7 \
  --dynamic-objective critical-gap-tail-integral \
  --minimum-dyadic-gap 3 \
  --dynamic-optimizer gradient --gradient-method lbfgs \
  --sobolev-order 4 --sobolev-cap 100 --threads 12
```

Here dyadic gap `m` means `high/low > 2^m`; the selected objective sums every
gap at least `m`. Its first projection-controlled run is under
`proof/l4/analysis/gap-tail-min3-h4-cap100/`.

Long searches can resume from a saved state without repeating lower cutoffs:

```bash
./build/navier_stokes_lab adversary --cutoffs 6 \
  --dynamic-warm-state proof/l4/states/gradient-continuation/dynamic/K6.tsv \
  --dynamic-generations 32 --dynamic-optimizer gradient \
  --state-dir proof/l4/states/gradient-continuation-refined
```

`--dynamic-objective q-gain` instead maximizes
`log(Q(T)/Q(0))`. This removes the reward for merely starting with a large `Q`
and targets relative time amplification, but it is ill-conditioned when
`Q(0)` approaches zero. The nondegenerate default for growth experiments is
`q-increase = Q(T)-Q(0)`. `terminal-q` is also available when the endpoint
value itself is the desired objective.

The self-test validates each derivative layer independently. Current relative
errors on the deterministic test state are approximately `2.6e-19` for FFT
JVP versus direct JVP, `1.9e-19` for FFT VJP versus direct VJP, `5.6e-19` for
RK4 adjoint duality, `2.0e-13` for a three-step `Q` trajectory gradient, and
`1.1e-12` for the critical trapezoidal time-integral gradient versus central
differences. A short `K=5` smoke certificate with 1,330 modes is stored in
`proof/l4/adversary/l4s-fft-gradient-K5-smoke.json`.

The first direct L4-A adjoint continuation is replayable with:

```bash
./build/navier_stokes_lab adversary \
  --cutoffs 1,2,3,4,5,6,7,8 \
  --dynamic-generations 12 \
  --dynamic-objective critical-integral \
  --dynamic-optimizer gradient \
  --sobolev-order 4 --sobolev-cap 100 \
  --nu 0.1 --evolve-time 0.01 --dt 0.002 \
  --threads 12 --backend auto \
  --state-dir proof/l4/states/critical-integral-h4-cap100
```

For this short-time, unit-energy experiment the optimized integral approaches
`6.01e-6` by `K=8`; a continued `K=8` run reaches `6.045e-6`. The final state
has squared homogeneous `H4` norm `15.76`, only modes through shell three are
active at relative energy `1e-8`, and the cutoff-shell energy is `1.49e-11`.
This is finite evidence for a smooth cutoff limit, not a uniform-in-time bound.
The certificate and interpretation are under
`proof/l4/analysis/critical-integral-h4-cap100/`.

A later same-state control found that 16 projected L-BFGS iterations improved
the integral about `34.5x` more than 16 steepest-gradient iterations. Projecting
the resulting `K=8` state to `K=3`, reoptimizing, and lifting it through
`K=3,...,8` gives refined integrals from `6.11840e-6` to `6.11975e-6`.
The increments after `K=4` are at most `2.0e-11`, while the top-shell energy
falls to `1.4e-13` at `K=8`. This is a reproducible smooth finite-cutoff branch,
not an analytical uniform bound.

## Fixed smooth projective family

Uniform-in-cutoff candidates should be tested on projections of one fixed
smooth initial field, not on unrelated random fields. The `family` command uses
deterministic per-wave coefficients with exponential Fourier decay and a common
normalization across all cutoffs:

```bash
./build/navier_stokes_lab family \
  --cutoffs 1,2,3,4 --decay 0.8 --nu 0.1 --time 0.1 --dt 0.002 \
  --seed-count 12 \
  --threads 12 --backend auto \
  --certificate proof/l4/family/l4s-family-scan.json
```

It verifies exact nesting, repeats every trajectory with half the time step,
and reports convergence of the critical integral, local/nonlocal stretching,
vorticity maximum, and high-vorticity direction coherence.

Both proof-oriented runners accept `--backend direct|fft|auto`. The built-in
radix-2 FFT uses a dealiased grid and is checked against direct triad summation
in every self-test. Optimizer restarts and small independent cutoffs run through
the CPU worker pool. A large single-cutoff FFT trajectory instead parallelizes
its component transforms and physical-grid diagnostics, avoiding nested worker
pools. The default limit is 12 CPU workers; use `--threads 1` for deterministic
performance debugging.

## One supporting simulation

Taylor–Green initial data:

```bash
./build/navier_stokes_lab simulate \
  --n 24 --nu 0.01 --time 1.0 --max-dt 0.002 \
  --init taylor-green --csv run.csv
```

Random Fourier modes, projected to the discrete divergence-free subspace:

```bash
./build/navier_stokes_lab simulate \
  --n 24 --init random --seed 7 --modes 12 --time 0.5
```

Diagnostics include kinetic energy, enstrophy, velocity and vorticity maximum
norms, divergence error, a Prodi–Serrin proxy, a BKM proxy, and the constant in
the tested discrete energy lemma.

## Resolution-based falsification

```bash
./build/navier_stokes_lab mine \
  --grids 12,16,24,32 --seeds 1,2,3,4 \
  --nu 0.01 --time 0.5 --report lemma-report.md
```

`mine` repeats runs for every `(grid, seed)` pair and measures whether peak
enstrophy, maximum vorticity, and the Prodi–Serrin integral grow with the
cutoff. A `candidate survives` result only means that this finite test set did
not falsify it.

## Discrete energy lemma

The solver uses skew-symmetric advection and an L2-orthogonal pressure
projection. One explicit Euler step therefore satisfies

```text
E[n+1] - E[n] + ν Δt ||∇u[n]||₂²
  <= (Δt² / 2) ||-N_h(u[n]) + ν Δ_h u[n]||₂².
```

The normalized theoretical constant is at most `1`; every time step checks it.
This is a stability lemma for the finite-dimensional scheme, not a global
regularity lemma for smooth 3D Navier–Stokes solutions.

## Tests

```bash
./build/navier_stokes_lab self-test
```

The tests cover the Helmholtz projection, divergence reduction, time stepping,
the discrete energy inequality, exact rational scaling, spectral triads,
direct/FFT agreement, JVP/VJP duality, central-difference gradient checks,
checkpointed reverse RK4, the shared local/nonlocal triad classifier, and
fixed-energy projected L-BFGS ascent with a steepest-gradient recovery path.
