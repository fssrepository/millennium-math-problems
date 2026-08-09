# Current Navier–Stokes research handoff

Updated: 2026-08-09. This is the post-audit resume point for the Codex
campaign archived as session `019fbc49-c7f3-7051-8627-e920169afd1b`.

## Resume here

1. Read the post-audit status at the top of
   [`proof/l4/lemmas/shifted-local-density/PALINSTROPHY_NORMALIZATION_TAIL.md`](proof/l4/lemmas/shifted-local-density/PALINSTROPHY_NORMALIZATION_TAIL.md),
   then the current-state section in [`PROOF_PLAN.md`](PROOF_PLAN.md).
2. Do **not** resume unrestricted PNT-2/4/5/7/8/12/14 optimization or treat
   PNT-12 as an active lemma. Under `u -> alpha u`, the open PNT quotients
   have amplitude degree `-1` and the squared Cauchy/Schur quotients have
   degree `-2`. The repository contains states with nonzero numerators, so
   their unrestricted suprema diverge as `alpha -> 0`.
3. The stored PNT searches remain valid evaluations of fixed-energy
   objectives. They do not test the withdrawn unrestricted statements because
   energy normalization removes the amplitude ray from the admissible set.
4. Before further PNT computation, re-derive the complete PNT-4/PNT-5 chain
   with an explicit admissible set or a justified homogeneous normalization.
   If no such formulation still closes the local/transition estimate, record
   the branch as failed and generate a different lemma candidate.
5. Independently, replace the floating-point sign check for the positive
   instantaneous local-density derivative by exact or interval validation.

Fast scope regression after building `navier_stokes_lab`:

```bash
./build/navier_stokes_lab self-test
```

The output must contain a passing `PNT amplitude/scope test` with open
quotient degree `-1`, Schur quotient degree `-2`, and a fixed-energy search
that does not cover the amplitude ray.

## Scientific checkpoint

- The Clay problem is **not solved**. L4 remains open.
- The moving-gap dynamic far-tail component has a cutoff-independent
  repository derivation. It is supplied for expert review; the gap-zero local
  block, logarithmic transition band, and full L4 lemma remain open.
- The unrestricted PNT chain as printed is withdrawn by exact amplitude
  scaling. No corrected PNT lemma is currently claimed.
- Fixed-energy PNT stress values, including the H16/K12 value
  `4.18215e-4`, remain reproducible finite-dimensional objective values. They
  are not evidence for an unrestricted analytical supremum.
- PNT-13 standalone height-shell decorrelation has strong finite-cutoff
  counterevidence. A uniform nonexistence theorem would still require a
  scalable counterexample family.
- Stored K3--K6 states give a positive directly evaluated instantaneous
  derivative of the local critical density in floating-point arithmetic.
  This is a numerical counterexample candidate until exact or interval
  arithmetic encloses the sign.
- The group-index `vjp_sum` path reduced K12 peak RSS from `7.94 GiB` to
  `4.88 GiB` without changing the fixed-energy objective and preserves
  replayable deterministic traces.

## Decision rule for the next attempt

- Put the admissible set next to every supremum and verify that the optimizer
  enforces the same set.
- Run the exact amplitude-degree audit before any numerical campaign. More
  fixed-energy compute cannot test a missing unrestricted amplitude quantifier.
- Do not insert an energy factor merely to repair dimensions. Re-derive any
  replacement from the PNT-4/PNT-5 closure chain and verify that it is still
  sufficient for the local/transition problem.
- If a homogeneous replacement survives the symbolic gate, then resume
  cutoff/height adversarial stress on precisely that objective.
- Append a failed-candidate row only with a clear status: exact analytical
  obstruction, rigorously validated finite-dimensional counterexample, or
  numerical counterevidence pending validation.
- Exact-gradient means the exact derivative of the implemented floating-point
  Galerkin/RK4 objective. It does not mean exact arithmetic or a result for the
  continuous PDE.

## Provenance and reading boundary

The private archive is
`~/.codex/archived_sessions/rollout-2026-08-01T09-46-19-019fbc49-c7f3-7051-8627-e920169afd1b.jsonl`
(about 124 MB; message span 2026-08-01 through 2026-08-03). It is not copied
into this repository. Open it only for an explicit provenance or timestamp
audit. The research state is carried by this handoff, `PROOF_PLAN.md`,
`proof/README.md`, the PNT note, and the candidate-status ledger.
