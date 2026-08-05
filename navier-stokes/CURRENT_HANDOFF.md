# Current Navier–Stokes research handoff

Updated: 2026-08-05. This is the short resume point for the Codex campaign
archived as session `019fbc49-c7f3-7051-8627-e920169afd1b`.

## Resume here

1. Read the current-state section in [`PROOF_PLAN.md`](PROOF_PLAN.md), then the
   stress ledger and reproduction block in
   [`proof/l4/lemmas/shifted-local-density/PALINSTROPHY_NORMALIZATION_TAIL.md`](proof/l4/lemmas/shifted-local-density/PALINSTROPHY_NORMALIZATION_TAIL.md).
2. Resume the **joint PNT-12 palinstrophy-normalization bound**, coupling the
   projective diagonal count and Gram row to the falling
   `S_full^2 B1/(Z^3 P^5)` normalization. Do not return to standalone shell
   decorrelation as the main route.
3. Computationally, extend/check the row-wise H32/H64 height/cutoff stress at
   K8/K12 before spending on a broader scan. Analytically, seek a
   cutoff-independent estimate for the coupled quantity, not separate loose
   bounds on its factors.

Fast certificate scan after building `navier_stokes_lab`:

```bash
./build/navier_stokes_lab local-sld-projective-normalization-cauchy-scan \
  --height 8 --state proof/l4/states/local-projective-normalization-cauchy/H8-K12-one-step/K12.tsv \
  --height 16 --state proof/l4/states/local-projective-normalization-cauchy/H16-K12-one-step/K12.tsv \
  --height 32 --state proof/l4/states/local-projective-normalization-cauchy/H32-K12-one-step/K12.tsv \
  --height 64 --state proof/l4/states/local-projective-normalization-cauchy/H64-K12-one-step/K12.tsv \
  --selection double-triple-remainder-without-123 --threads 12 \
  --certificate /tmp/K12-normalization-cauchy-height-scan.json
```

## Scientific checkpoint

- The Clay problem is **not solved**. L4 remains open.
- `F000`–`F015` are sixteen rejected lemma candidates. The canonical ledger is
  [`proof/failed_lemmas.tsv`](proof/failed_lemmas.tsv); do not rediscover them.
- One dynamic far-tail component is proved with cutoff-independent constants.
  `L4.1b`, the local/transition closure and the full L4 lemma remain open; see
  [`proof/README.md`](proof/README.md).
- PNT-13 standalone height-shell decorrelation has strong finite-cutoff
  counterevidence: raw correlations reach roughly `0.9857–0.9990` and weighted
  constants track the cutoff-wall maximum. A uniform nonexistence claim still
  needs a scalable counterexample family, so record it as a failed mechanism,
  not a theorem that every related joint bound fails.
- PNT-12 survives as the active coupled target. Direct exact-gradient search
  raised the finite H16/K12 record from `3.89149e-4` to `4.18215e-4` (`+7.46%`).
  That disproves the old apparent decimal plateau, not the inequality.
- The group-index `vjp_sum` path reduced K12 peak RSS from `7.94 GiB` to
  `4.88 GiB` without changing the objective and preserves replayable
  deterministic traces.

## Decision rule for the next attempt

- If height/cutoff continuation produces persistent growth, isolate its shell
  geometry and try to turn it into a projective scalable counterfamily.
- If the coupled values stabilize, stop extrapolating decimals and prove a
  cutoff-independent joint lemma retaining diagonal, Gram-row and normalization
  tradeoffs.
- Append a new `Fxxx` row only after a concrete obstruction is replayable.
- Exact-gradient means the exact derivative of the implemented finite
  Galerkin/RK4 objective. It does not mean an exact result for the continuous
  PDE.

## Provenance and reading boundary

The private archive is
`~/.codex/archived_sessions/rollout-2026-08-01T09-46-19-019fbc49-c7f3-7051-8627-e920169afd1b.jsonl`
(about 124 MB; message span 2026-08-01 through 2026-08-03). It is not copied
into this repository. Open it only to audit exact wording/timestamps. The
research state is carried by this handoff plus `PROOF_PLAN.md`, `proof/README.md`,
the lemma note above and the certificate/failed-lemma ledgers.
