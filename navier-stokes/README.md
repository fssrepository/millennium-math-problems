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

- `src/app/` contains only command dispatch and the integrated proof runner;
- `src/spectral/` contains Fourier states, direct/FFT Galerkin dynamics,
  objectives, adjoints, and trajectory evaluation;
- `src/optimization/`, `src/triads/`, `src/helical/`, and `src/proof/` isolate
  their respective algorithms;
- `src/local_sld/core/`, `analysis/`, `optimization/`, and `cli/` separate the
  active shifted-local-density lemma work into numerical primitives,
  diagnostic ledgers, searches, and artifact/report handling;
- `src/reporting/` contains the remaining shared certificate writers.

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
  exact `Q=D^4 Z`, `integral D^4 Z^2 dt`, endpoint critical-density increase,
  and shifted log-gain gradients, checkpointed reverse RK4 passes, and
  fixed-energy Riemannian gradient search;
- `DynamicAdversary` owns one complete forward/adjoint optimization context,
  while `DynamicAdversaryEnsemble` runs independent contexts concurrently and
  refines only the winning trajectory with a halved time step;
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
- `HelicalGapLedger` resolves those sectors by dyadic gap, while
  `LocalTriadSymmetrizer` groups complete signed Fourier triads and certifies
  the exact frequency-spread bound coming from triadwise energy cancellation;
- `OrthogonalTriadGeometry` proves and enumerates the bounded-degree lattice
  graph for equal-length orthogonal waves and certifies its subcritical
  `K^(-1/2)` transfer-to-viscosity scaling;
- `LocalSignatureGeometry` extends the degree bound to every fixed
  squared-length triple, square-sums the signature transfers, and certifies
  that an effective-signature exponent below one is viscosity-subcritical;
- `LocalSignatureObjective` supplies central-difference-certified analytic
  gradients for signed amplification and absolute local transfer;
- `LocalSignatureDensity` owns the exact coupled critical factorization, while
  `LocalSignatureFactorAdversary` falsifies pointwise factor-correlation
  mechanisms with parallel one-step RK4 probes;
- `ShiftedCriticalDensityLemma` certifies the exact homogeneity and
  conditional Gronwall closure of the scale-compatible `E(0)P(0)` shift;
  `ShiftedCriticalDensityAnalyzer` evaluates its exact instantaneous
  Navier--Stokes derivative without an RK4 horizon;
- `LocalCriticalDerivativeLedger`, `StretchingDerivativeLedger`, and
  `ShiftedCriticalDensityBudgetAnalyzer` split that derivative by `S/Z/P`
  chain-rule source, Euler/viscous dynamics, local/nonlocal RHS, and all three
  Frechet slots of `<A u,B_local(u,u)>`, with independent reconstruction and
  PDE-identity residuals;
- `LocalQuarticIdentityLedger` preserves the exact negative-square outer-slot
  identity and the local/nonlocal enstrophy-transfer identities before any
  inequality discards their signs;
- `LocalQuarticShellLedger` resolves the complete normalized local quartet by
  target mode, exact Laplace eigenshell, hard shell, and dyadic annulus;
- `LocalQuarticShellEnvelope` proves and certifies the cutoff-independent
  bound `||A^(1/2)B_local||_2^2 <= 7264120.5 sqrt(Z) P^(3/2)` using only
  elementary lattice counting, shell overlap, and `H^(3/2)` interpolation;
- `LocalQuarticCommutator` combines the outer and advected Frechet slots into
  the exact multiplier `|k|^2-|q|^2` and checks its duality and symbol bounds;
- `LocalQuarticProjectedResidual` combines the outer, enstrophy, and
  palinstrophy chain-rule entries and verifies their exact completed-square
  representation before the remaining transport slots are estimated;
- `LocalQuarticReducedLedger` merges both exact reductions and emits the
  two-entry, denominator-free local SLD-1P numerator with an independent
  polynomial reconstruction check;
- `LocalQuarticClosureTarget` certifies the exact weighted-Young implication
  from the explicit `K+G` depletion target to the local SLD-1P inequality;
- `LocalQuarticClosureObjective` reverse-differentiates both the sufficient
  absolute closure ratio and the actual signed local SLD-1P quotient;
- `LocalQuarticClosureAdversary`, its CLI, and its reporter run 12-worker
  projected L-BFGS cutoff continuations and save English JSON plus replayable
  Fourier states without placing artifacts at an analysis/adversary root;
- `LocalSldCyclicBasis` owns the reusable cyclic axis-shear and quadratic
  response basis, while `LocalSldCyclicAnsatz` optimizes its static mixture;
- `LocalSldCyclicTrajectoryAnsatz` searches the frozen-data trajectory
  maximum in that one-angle family, refines the time step, and uses the exact
  adjoint to measure both the restricted and full projected gradients;
- `LocalSldCyclicKrylovAnsatz` adds the orthogonalized advection-JVP response
  and optimizes its three coefficients on the fixed-energy sphere;
- `LocalSldResponseHierarchy` constructs the full quadratic response
  recursion, while `LocalSldCyclicOrbitBasis` supplies the missing transverse
  `(2,1,1)` polarization and both oriented `(3,1,0)` cyclic orbits;
- `LocalSldResponseDiagonal` excludes Galerkin-wall-contaminated response
  orders, while `LocalSldResponseTensor` evaluates every exact direct-triad
  interaction and its orthogonal complement with separate input/output
  analytic radii;
- `LocalSldResponseBasis` performs degree-ordered orthonormalization of scalar
  responses and transverse cyclic orbits; the tensor can greedily insert
  missing quadratic products without contaminating low analytic degrees with
  higher-cutoff response directions;
- `LocalSldTwoScaleState` constructs exact cyclic low/high dilations without
  repeated response generation; `LocalSldDoublingScaleScan` evaluates the
  doubling, remainder, and mixed blocks across scale and response angle;
- `LocalSldDoublingShellLedger` emits the exact ordered dyadic-shell matrix,
  while `LocalSldProjectedSquare` verifies the completed-square identity for
  any fixed triad selection;
- `LocalSldTrajectoryEvaluator` evaluates and dt-refines any saved state with
  either the direct RK4 oracle or the FFT forward/VJP backend;
- `LocalSldSignatureBlock` splits `K+G` exactly into a selected squared-length
  signature family, its local complement, and their independently evaluated
  mixed terms;
- `LocalSldBlockObjective` gives those three terms exact gradients under the
  common full-SLD normalization, so their real contributions can be searched
  without replacing the full stretching by a block stretching;
- `LocalSldTrajectoryAdjoint` freezes `k0,B0` from the initial datum, evolves
  the Galerkin state, and reverse-differentiates the terminal local quotient
  through RK4 and through both frozen parameters;
- `LocalSignatureGradientAdversary` runs 12-worker targeted counterexample
  searches, while `LocalSignatureTrajectoryAnalyzer` verifies the exact
  dynamic signature factorization across Galerkin cutoffs;
- `HelicalSectorObjective`, `HelicalSectorAdjoint`, and the two helical
  adversaries provide exact static and checkpointed trajectory gradients for
  sector-selective local searches;
- `HelicalAdversaryCli` runs replayable fixed-energy sector searches with
  parallel multistarts, while `HelicalCutoffScan` evaluates one identical
  initial state on several Galerkin cutoffs and emits an aggregate convergence
  certificate;
- `MovingGapController` implements the logarithmic gap choice that turns the
  post-Young far-tail remainder from `Z^3` into a linear `Z` term;
- `DyadicShellBounds` isolates the cutoff-independent scalar shell sums,
  including explicit `sqrt(2)`, `sqrt(8/7)`, and Holder constants;
- `PeriodicShellGeometry` and `PeriodicTailBound` supply explicit torus
  lattice, adjacent-shell, and final hard-shell FT-1 constants;
- `FarTailClosure` verifies the cutoff-independent moving-gap Young reduction,
  while `TransitionBlockScaling` rejects logarithmic band counting as a
  closure mechanism;
- `DoublingQuartetClosure` combines the equal-length orthogonal incidence
  bound with exact rational quartet power counting: the complete doubling
  family now closes at the cutoff-independent `Z^(5/4)P^(3/4)` scale; an exact
  two-scale counterexample is retained because it rejects a tempting but
  unnecessary intermediate absorption inequality;
- `RemainderQuartetClosure` computes the exact dense-signature loss and the
  effective `R^(3/2)` incidence degree required by any unsigned route;
- `LocalSldRemainderSignatureLedger` reconstructs the closed remainder
  bracket signature by signature in one parallel interaction pass;
- `LocalSldRemainderDoubleSquare` removes the raw dense advection norm from
  the one-sided target by two exact square completions, while
  `LocalSldRemainderEnvelopeObjective` supplies the exact VJP for the
  resulting commutator envelope;
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
report generation are separate compilation units. Pointwise bounds for both
raw signature participation and signed amplification have been adversarially
rejected. The dominant doubling quartet now has a conventional
cutoff-independent proof; the next mathematical tasks are the closed local
signature remainder and its mixed interaction with the doubling block. The
full local SLD lemma and the Clay problem remain open. This keeps rebuilds
dependency-free and makes each layer independently replaceable.

The current direct local-lemma search is reproducible with:

```bash
./build/navier_stokes_lab doubling-quartet-certificate \
  --max-cutoff 12 \
  --certificate proof/l4/analysis/shifted-local-density/doubling-quartet/closed-family-K12.json

./build/navier_stokes_lab local-sld-doubling-scale-scan \
  --min-scale 2 --max-scale 12 \
  --angle-min -1.2 --angle-max 1.2 --angle-count 25 \
  --energy-decay-power 2.75 --threads 12 \
  --certificate proof/l4/analysis/shifted-local-density/doubling-quartet/two-scale-L2-L12-angle-scan.json

./build/navier_stokes_lab local-closure-adversary \
  --objective sld-ratio --min-cutoff 1 --max-cutoff 4 \
  --restarts 12 --workers 12 --iterations 20 --method lbfgs \
  --certificate proof/l4/adversary/shifted-local-density/direct-local-sld-ratio-K1-K4.json \
  --state-dir proof/l4/states/local-sld-ratio/K1-K4

./build/navier_stokes_lab local-sld-ansatz \
  --samples 8192 --refinements 112 \
  --certificate proof/l4/analysis/shifted-local-density/cyclic-two-basis-ansatz.json \
  --state proof/l4/states/local-sld-ratio/cyclic-ansatz/K2.tsv

./build/navier_stokes_lab local-sld-trajectory-ansatz \
  --cutoff 2 --samples 256 --refinements 64 --threads 12 \
  --trajectory-steps 500 --nu 0.1 --dt 0.001 --backend auto \
  --certificate proof/l4/analysis/shifted-local-density/cyclic-trajectory-T050-K2.json \
  --state proof/l4/states/local-sld-trajectory/cyclic-ansatz-T050-K2/K2.tsv

./build/navier_stokes_lab local-sld-response-hierarchy \
  --state proof/l4/states/local-sld-trajectory/maximum-T050-K3-from-K2/K3.tsv \
  --depth 16 --threads 12 --include-211-transverse --include-310-orbits \
  --projected-state proof/l4/states/local-sld-trajectory/response-hierarchy-projection-K3/depth16-plus-orbits.tsv \
  --residual-state proof/l4/states/local-sld-trajectory/response-hierarchy-residual-K3/depth16-plus-orbits.tsv \
  --certificate proof/l4/analysis/shifted-local-density/cyclic-response-hierarchy-K3-depth16-plus-orbits.json

./build/navier_stokes_lab local-sld-response-diagonal \
  --state proof/l4/states/local-sld-trajectory/maximum-T050-K2/K2.tsv \
  --state proof/l4/states/local-sld-trajectory/maximum-T050-K3-from-K2/K3.tsv \
  --state proof/l4/states/local-sld-trajectory/maximum-T032-K4-from-K3/K4.tsv \
  --max-depth 16 --radius 1.25 --threads 12 \
  --certificate proof/l4/analysis/shifted-local-density/cyclic-response-diagonal-K2-K4-r125.json

./build/navier_stokes_lab local-sld-response-tensor \
  --cutoff 12 --depth 13 \
  --input-radius 2 --output-radius 1.15 \
  --tolerance 1e-14 --threads 12 \
  --certificate proof/l4/analysis/shifted-local-density/response-tensor/R200-r115/K12.json

./build/navier_stokes_lab local-sld-response-tensor \
  --cutoff 5 --depth 6 --input-radius 2 --output-radius 1.15 \
  --include-211-transverse --include-310-orbits \
  --closure-extensions 16 --tolerance 1e-14 --threads 12 \
  --certificate proof/l4/analysis/shifted-local-density/response-tensor/augmented-graded-R200-r115/K5-closure16.json

./build/navier_stokes_lab local-sld-trajectory-evaluate \
  --state proof/l4/states/local-sld-trajectory/response-hierarchy-projection-K3/depth16-plus-orbits.tsv \
  --trajectory-steps 320 --threads 12 --nu 0.1 --dt 0.001 --backend auto \
  --certificate proof/l4/analysis/shifted-local-density/cyclic-response-hierarchy-K3-depth16-plus-orbits-trajectory.json

./build/navier_stokes_lab local-sld-block \
  --state proof/l4/states/local-sld-ratio/K8/K8.tsv \
  --doubling-family --threads 12 \
  --certificate proof/l4/analysis/shifted-local-density/direct-local-sld-K8-doubling-family.json

./build/navier_stokes_lab local-closure-adversary \
  --objective block-ratio --selection doubling-family \
  --min-cutoff 1 --max-cutoff 6 --restarts 12 --workers 12 \
  --iterations 30 --method lbfgs \
  --certificate proof/l4/adversary/shifted-local-density/common-sld-doubling-block-K1-K6.json \
  --state-dir proof/l4/states/local-sld-block/doubling-family-K1-K6

./build/navier_stokes_lab local-closure-adversary \
  --objective terminal-sld-ratio --selection local \
  --min-cutoff 1 --max-cutoff 3 --restarts 12 --workers 12 \
  --iterations 20 --trajectory-steps 10 --nu 0.1 --dt 0.001 \
  --certificate proof/l4/adversary/shifted-local-density/frozen-terminal-sld-T001-K1-K3.json \
  --state-dir proof/l4/states/local-sld-trajectory/frozen-T001-K1-K3

./build/navier_stokes_lab local-closure-adversary \
  --objective maximum-sld-ratio --selection local \
  --min-cutoff 3 --max-cutoff 3 --restarts 12 --workers 12 \
  --iterations 8 --trajectory-steps 500 --nu 0.1 --dt 0.001 \
  --backend auto \
  --warm-state proof/l4/states/local-sld-trajectory/maximum-T050-K2/K2.tsv \
  --certificate proof/l4/adversary/shifted-local-density/frozen-maximum-sld-T050-K3-from-K2.json \
  --state-dir proof/l4/states/local-sld-trajectory/maximum-T050-K3-from-K2
```

Trajectory searches accept `--backend direct`, `--backend fft`, or
`--backend auto`. `direct` remains the reference oracle; `auto` switches the
forward RK4 and its exact discrete VJP to the validated FFT path from K3 while
the selected local closure objective retains exact triad sums.
Each completed optimizer restart is also written immediately below
`STATE_DIR/restarts/K*/R*.tsv`. A machine or process interruption can therefore
lose only the currently running restarts, not every completed branch.

The first command maximizes the signed polynomial quotient itself. The
`closure-ratio` objective remains available as a stronger sufficient screen,
but it is not used as a substitute for the direct target.
The common-normalization scans stabilize near `7.72930e-4` for the doubling
family, `2.20683e-4` for its complement, and `1.34936e-4` for the mixed term.
A separate signed-bracket search rejects the shortcut `K+G<=0` with the
positive value `0.0740187069851`; that counterexample has nearly zero
stretching, which is why the code keeps the full joint SLD factor.
The canonical cyclic axis shear exactly attains the static absolute value
`|c|=1/3`, and padded 12-start searches find no larger value through K8.
This is a sharp conjecture, not a proof. With `k0,B0` frozen, K3 terminal
optimization increases from `7.98918e-4` at `T=0.01` to `8.46863e-4` at
`T=0.20`; this is now the active trajectory-level target.
The maximum-on-trajectory version removes endpoint bias. Warm continuation
from the K2 winner corrects an earlier secondary branch and gives the refined
K3 lower bound `8.53498799310e-4` at `t=0.298` inside `[0,0.5]`; dt halving
changes it by `8.16e-16` relatively. At that checkpoint, the exact frozen
block split is `8.24077287430e-4 + 1.51420260355e-5 + 1.42794858448e-5`,
with `6.20e-20` reconstruction error. The doubling family supplies `96.55%`
of the absolute three-block sum. The winning state retains `99.969%` of its
energy in the first hard shell, but its projected gradient is still
`4.04e-4`; this remains a lower bound, not a converged global maximum.

Six accepted K4 continuation steps improve the finite lower bound to
`8.53527437357e-4` at the same `t=0.298`; the dt-halving error is `9.00e-16`.
The added shell-four energy is `1.97e-9`, the K4-to-K3 state residual is
`2.96e-4`, and the final projected objective gradient is `5.17e-5`.
At the refined peak the exact three-block split is
`8.24057223723e-4 + 1.51705418626e-5 + 1.42996717717e-5`.

The explicit response hierarchy makes that pattern quantitative. Orders
`0..15`, the transverse `(2,1,1)` polarization, and both oriented `(3,1,0)`
orbits capture `99.9998628%` of the K3 winner's energy. Evolving only its
19-vector projection gives `8.53172846573e-4` at the same `t=0.298`, which is
`99.96181%` of the unrestricted lower bound. Direct and FFT trajectory
evaluation agree to all reported digits. This identifies a compact candidate
extremal structure; it is not a cutoff-uniform theorem.

On the K4 winner the same 19 directions capture `99.9997617%` of state energy
and `99.9331554%` of the trajectory objective. The cutoff-diagonal response
ledger then retains only the boundary-free orders `0..K`. On the K4 state
zero-padded through K8, `A_1.25(K)` stays below `1.15916`; the full radius
sweep and the exact sequence majorant are recorded in
`proof/l4/lemmas/shifted-local-density/RESPONSE_DIAGONAL.md`. This is the
current route toward a weighted response-space lemma, not a proof of its
required cutoff-uniform operator or complement bounds. The exact interaction
tensor sharpens the next target: with input radius `2` and output radius
`23/20`, every ordered pair through K12 satisfies
`sum_m (23/20)^m |<b_m,B(b_i,b_j)>| <= [23/(20 sqrt(3))]2^(i+j)`, with
equality at the axis pair. Proving it for all response orders, together with
a shell-resolved transverse-complement estimate, is the current analytical
task. In the degree-graded augmented tree, sixteen closure extensions give
projected block constants `1.04467`, `1.04894`, and `1.04894` at K3--K5; the
corresponding projected-plus-complement finite bounds are `1.17088`, `1.30595`,
and `1.32010`. See
`proof/l4/lemmas/shifted-local-density/RESPONSE_TENSOR.md`.

For the universal dominant block, equal-length orthogonal incidence proves the
one-shell bound `|K_d+G_d|_j <= C R_j^5 E_{near,j}^2`, which is half a
derivative better than the LQC-3 target. Neighbor-shell locality sums the
structural entries. Direct sequence estimates give
`S <= C Z^(5/4)P^(1/4)` and `T <= C Z^(1/4)P^(5/4)`, closing both global
normalization terms at `Z^(5/4)P^(3/4)`. Thus the complete doubling block is
now proved cutoff-independently. The active proof targets are the closed
signature remainder and the mixed block. For the remainder, the signed LQC-3
search is flat at `0.102959` through K8, whereas a growing absolute branch is
negative and is generated by the coercive square. Two exact square
completions reduce the real open statement to a commutator/negative-square
absorption inequality; bounding the positive commutator envelope alone is
too strong. See
`proof/l4/lemmas/shifted-local-density/DOUBLING_QUARTET.md` and
`proof/l4/lemmas/shifted-local-density/REMAINDER_QUARTET.md`.

The helical local target has its own replayable optimizer and same-state
cutoff scan:

```bash
./build/navier_stokes_lab helical-adversary \
  --state proof/l4/states/helical-heterochiral-local-integral/dynamic/K3.tsv \
  --state-output /tmp/helical-K3.tsv --certificate /tmp/helical-K3.json \
  --selection heterochiral --spread broad --mode trajectory \
  --cutoff 3 --iterations 8 --trajectory-steps 10 \
  --restarts 12 --workers 12 --restart-mutation 0.03 \
  --nu 0.1 --dt 0.001

./build/navier_stokes_lab helical-cutoff-scan \
  --state /tmp/helical-K3.tsv --certificate /tmp/helical-scan.json \
  --selection heterochiral --min-cutoff 3 --max-cutoff 6 \
  --trajectory-steps 10 --workers 12 --nu 0.1 --dt 0.001
```

Setting `--iterations 0` evaluates an adversary state without optimizing it,
which separates cutoff effects from optimizer progress.
The spread mask is `all`, `equal`, `narrow` (relative squared-frequency spread
at most `1/4`), or `broad` (its gap-zero complement).

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
  --dynamic-restarts 1 \
  --dynamic-generations 24 --dynamic-objective max-q \
  --dynamic-optimizer gradient \
  --nu 0.1 --evolve-time 0.1 --dt 0.002 \
  --certificate proof/l4/adversary/gradient/l4s-maxq-adversary.json \
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

`--dynamic-restarts N` launches genuinely independent dynamic adjoint
optimizations. Additional starts alternate between Sobolev-retracted
perturbations of the warm continuation and retracted random states. Their
final objectives and the winning restart are recorded in the JSON artifact;
the winning state alone receives the refined `dt/2` trajectory check.

The two terms of the local/nonlocal proof split can be attacked directly:

```bash
./build/navier_stokes_lab adversary \
  --cutoffs 3 --dynamic-generations 16 \
  --dynamic-objective critical-nonlocal-integral \
  --dynamic-optimizer gradient --gradient-method lbfgs \
  --sobolev-order 4 --sobolev-cap 100 \
  --nu 0.1 --evolve-time 0.01 --dt 0.002 \
  --threads 12 \
  --certificate proof/l4/adversary/partitioned-integral/nonlocal.json
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
`proof/l4/adversary/gradient/l4s-fft-gradient-K5-smoke.json`.

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
