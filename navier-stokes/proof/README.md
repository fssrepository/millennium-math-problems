# Proof artifacts

Machine-generated certificates and replayable states are grouped by proof
stage:

- `l2/`: exact scaling and Fourier-triad obstruction certificates;
- `l4/adversary/`: static and trajectory adversary certificates;
- `l4/family/`: fixed smooth projective-family certificates;
- `l4/states/<experiment>/static/`: replayable static-search winners;
- `l4/states/<experiment>/dynamic/`: replayable trajectory-search winners;
- `l4/analysis/`: shell, Sobolev, active-mode, and projectivity certificates;
- `failed_lemmas.tsv`: the global rejection ledger.

JSON and TSV files are computational evidence or exact finite algebra checks.
Their presence does not replace the infinite-dimensional estimates required by
`PROOF_PLAN.md`.

`l4/adversary/l4s-fft-gradient-K5-smoke.json` is the first certificate produced
through the complete FFT forward/adjoint path. It records a fixed-energy
`max-q` search at cutoff `K=5` with 1,330 Fourier modes and is intended as a
fast regression artifact, not as evidence of asymptotic cutoff behavior.

`l4/adversary/l4-critical-integral-h4-cap100-K1-K8.json` is the first cutoff
continuation that directly maximizes the L4-A time integral with the exact
discrete adjoint. Its states are under
`l4/states/critical-integral-h4-cap100/`; the warm-start refinement is under
`l4/states/critical-integral-h4-cap100-refined/`. The numerical interpretation
and state-family certificate are in
`l4/analysis/critical-integral-h4-cap100/README.md`.

`l4/adversary/l4-critical-integral-h4-cap100-projective-lbfgs-K3-K8.json`
records the improved low-mode-to-high-mode projected L-BFGS continuation. Its
iteration traces distinguish optimizer progress from genuine cutoff gain. The
matching state-family analysis is under
`l4/analysis/critical-integral-h4-cap100-projective-lbfgs-K3-K8/`.

`l4/analysis/partitioned-critical-integrals/` records the first adversaries
that maximize the exact local and nonlocal L4-A pieces separately. It includes
the replayable `K=3` shell analyses and an initial-`H4` cap sweep. The data
invalidate the branch-specific assumption that the nonlocal term is always
negligible; they do not prove either partitioned lemma.

`l4/analysis/far-nonlocal-h4-cap100/` refines the nonlocal class into the first
transition band and the truly scale-separated tail. It contains the first
fine-time-step, projection-controlled `K=3,...,6` far-tail adversary.

`l4/analysis/gap-tail-min3-h4-cap100/` tests the next dyadic tail
`high/low > 8` with a parameterized exact objective. It records fine-time-step
checks and both downward projections and zero-padded upward cutoff controls.

`l4/lemmas/dynamic-far-tail/` contains the current human-checkable proof draft:
the three-role paraproduct envelope, the fixed-gap obstruction, and the
moving-gap Young closure that leaves the logarithmic transition/local block.

`l4/lemmas/transition-block/` records the exact scaling obstruction showing
that counting the remaining `O(log Z)` bands does not by itself close the
enstrophy inequality.

`l4/lemmas/local-helicity/` records the exact eight-sector helical triad
decomposition, the heterochiral dominance of current local extremizers, and
the pure-helicity counterexample to exact homochiral cancellation.
