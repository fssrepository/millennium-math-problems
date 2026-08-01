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
the pure-helicity counterexample to exact homochiral cancellation. The
replayable K3--K6 same-state continuation is summarized there and certified by
`l4/adversary/helical-heterochiral-cutoff-scan-K3-K6.json`.

`l4/lemmas/local-triad-symmetry/` proves the exact complete-triad
frequency-spread estimate `LS-1`. It eliminates equal-frequency enstrophy
transfer but records why broad-spread local triads still leave the key L4
depletion open.

`l4/lemmas/local-orthogonal-triads/` controls the infinite signature family
`(|p|^2,|q|^2,|p+q|^2)=(r,r,2r)`. Its bounded lattice degree improves the
local transfer from critical frequency power `9/2` to `7/2`, making the family
viscosity-absorbable at high frequency.

`l4/lemmas/local-signature-families/` proves the same `O(K)` lattice-degree
bound for every fixed squared-length signature and square-sums their
transfers. The full broad local block closes if its effective coherent
signature count is `O(K^mu)` for any `mu<1`; the replayed K3 broad endpoint has
`N_eff=1.1902068038`. The unrestricted pointwise hypothesis is rejected by
`l4/analysis/local-signature-adversary/`, whose flat-spectrum search reaches
`N_eff=1748.62` at K6. The sign-preserving amplification `A_sig` remains only
`2.33549` under random sampling, but the exact-gradient adversary in
`l4/analysis/local-signature-gradient/` raises it to `49.3891` and rejects the
pointwise hypothesis. The surviving target is the time integral of the exact
coupled factorization `A_sig^4 R^2/(Z P^3)`. A true 12-worker dynamic
multistart in `l4/analysis/local-signature-coupled-integral/` converges to the
same smooth K3--K6 branch and does not expose cutoff concentration.

`l4/analysis/local-critical-increase/` uses the exact discrete endpoint
gradient to reject monotonic decay of the coupled local density (F010). Its
positive-growth K3--K6 branch is time-step stable and spectrally convergent,
but its increase flattens with cutoff and therefore supplies a mechanism test,
not a blow-up certificate.

`l4/lemmas/shifted-local-density/` states the next explicit analytical
candidate. The initial shift `B0=E(0)P(0)` has exactly the same amplitude and
Navier--Stokes scaling as the local critical density. A cutoff-uniform bound on
`d log(C_local+B0)/dt` by `k0 Z` would close the local time integral through
the energy identity; this differential estimate remains open. The directory
also records the exact derivative ledger, its local/nonlocal and
Euler/viscous decomposition, and the denominator-free polynomial obligation
`SLD-1P` that is now the analytical restart point.
