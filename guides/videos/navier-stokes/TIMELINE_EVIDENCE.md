# Research timeline and screen mapping

## Scope

The main Navier–Stokes campaign began at `2026-08-01 07:51:01 UTC` and reached
the final reported PNT-12 stress result at `2026-08-02 17:19:22 UTC`:
**33 h 28 min 21 s**. It crossed two local calendar days in
Europe/Bratislava. The later session continued documenting and estimating the
open work, but it is not added to the research clock shown in the film.

The film maps the full 33:28:21 research interval onto `00:12.40–00:25.20`.
That 12.80-second window is linear: one screen second represents approximately
2 h 37 min. Event labels are rounded to the nearest minute.

## Audited milestones

| Real timestamp (UTC) | Elapsed | Film time | Evidence / event | Public shorthand |
| --- | ---: | ---: | --- | --- |
| 2026-08-01 07:51:01 | T+00:00 | 00:12.40 | Campaign begins | `START` |
| 2026-08-01 08:26:42 | T+00:35 | 00:12.63 | F002 rejected | `SCALING • REJECTED` |
| 2026-08-01 14:21:15 | T+06:30 | 00:14.89 | F006 homochiral-zero claim rejected | `CANCELLATION • REJECTED` |
| 2026-08-01 17:00:57 | T+09:10 | 00:15.91 | F009 stable simultaneous growth | `GROWTH • REJECTED` |
| 2026-08-01 20:25:32 | T+12:34 | 00:17.21 | F012 positive K+G counterexample | `POSITIVE CASE • REJECTED` |
| 2026-08-02 03:26:41 | T+19:36 | 00:19.89 | F013 rejected by exact-gradient search | `F013 • REJECTED` |
| 2026-08-02 04:37:33 | T+20:46 | 00:20.34 | F014 stronger dense branch | `F014 • REJECTED` |
| 2026-08-02 06:26 approx. | T+22:35 | 00:21.04 | Memory failure / saved-state resume | `CHECKPOINT → RESUME` |
| 2026-08-02 16:25:18 | T+32:34 | 00:24.86 | F015 / PNT-13 shell-orthogonality route rejected | `PNT-13 • REJECTED` |
| 2026-08-02 17:19:22 | T+33:28 | 00:25.20 | PNT-12 finite stress campaign reported; cutoff-independent claim remains open | `PNT-12 • OPEN FRONTIER` |

## Why the film inspects F015

F000–F015 are 16 different rejected candidates with different obstructions.
The short film does not claim that all of them failed because of one tail. It
uses F015 / PNT-13 as one representative failure that can be made legible on a
phone.

PNT-13 expected separated projective height-shell correlations to decay like
`2^(-gap)`. The exact-gradient adversaries instead reached raw correlations
`0.985710`, `0.998998` and `0.994971` at gaps 4, 6 and 7. The corresponding
weighted constants rose to `15.7714`, `63.9359` and `127.356`, close to the
respective `2^gap` maxima. The mobile overlay therefore uses the faithful
shorthand `FOUND • CORRELATION ≈ 1` and points at a persistent
**frequency-shell tail**. It is a concept visualization, not literal dye-tail
evidence from the C++ calculation.

That rejection removed standalone shell almost-orthogonality as the active
mechanism. The AI loop saved the counterexample and routed the next work to the
coupled PNT-12 quantity, which remains open.

## Final state shown in S07

- 16 rejected candidates: F000–F015.
- One partial far-tail estimate proved.
- The cutoff-independent L4 / PNT-12 step remains open.
- Finite stress record improved from `3.89149e-4` to `4.18215e-4`, about
  `+7.46%`; this belongs in a long description, not the mobile overlay.
- K12 memory was reduced from `7.94 GiB` to `4.88 GiB`.
- The Clay Millennium problem remains unsolved.

## Estimate wording

The later session estimate was:

- next PNT-12 decision: about 3–8 serious iterations, roughly 1–3 days of
  targeted calculation/analysis;
- if it survives larger H32/H64 and cutoff stress, an analytic
  cutoff-independent lemma may take days to weeks, with no guarantee;
- even a PNT-12 proof would leave the other channel, RQ-11 and L4→L5→L6 work;
- a complete Clay solution cannot be estimated honestly.

The short therefore says `~1–3 days` only for the **next decision**, not for a
proof or the full problem.

## Laptop and cluster interpretation

The observed campaign ran on a 13th-generation Intel Core i5-1335U laptop with
12 logical CPUs. The code already uses 12-worker restart pools, parallel
multistarts and parallel cutoff/interaction passes. A cluster can distribute
independent cutoffs, seeds, restarts and parameter sweeps across nodes.

For a concrete mobile comparison, the film uses a **10-node cluster** and says
`up to ~10×` for **parallel sweeps only**. This is an idealized engineering
ceiling inferred from ten comparable job groups running concurrently; it is
not a cluster benchmark. Practical speedup may be lower because of memory,
I/O, load imbalance and stages that remain serial. In particular:

- do not divide the 33 h 28 min campaign duration by ten;
- do not convert the `~1–3 days` next-decision estimate mechanically into
  hours;
- candidate formulation, AI interpretation, certificate review and analytic
  proof do not scale linearly with node count;
- a larger cluster changes search wall time, not the truth of a lemma or the
  absence of a full-solution estimate.

The exact on-screen qualifier is `COMPUTE ESTIMATE • TESTS ONLY`.

## Clay prize context

Clay established seven Millennium Prize Problems with US$1 million allocated
to each. The Poincaré Conjecture has been resolved, leaving six unresolved.
The reward is not automatic after an unrefuted forum post. Before CMI considers
a proposed solution, the official rules require:

1. publication in a Qualifying Outlet;
2. at least two years since that publication;
3. general acceptance in the global mathematics community.

CMI then decides whether the proposal merits its detailed consideration. The
film compresses this to `QUALIFYING PUBLICATION → ≥2 YEARS → GENERAL ACCEPTANCE
→ CMI REVIEW`.

## Provenance

Primary local record:

`~/.codex/archived_sessions/rollout-2026-08-01T09-46-19-019fbc49-c7f3-7051-8627-e920169afd1b.jsonl`

Rejected-lemma ledger:

`../../../navier-stokes/proof/failed_lemmas.tsv`

The private session archive is referenced for auditability but is not copied
into this production bundle.

## Authorship evidence and public wording

The archived campaign shows a connected AI research loop, not merely AI-written
boilerplate around a human-run solver. Under human direction, Codex repeatedly:

- turned the current proof gap into candidate statements and finite tests;
- wrote and refactored the C++ implementation;
- built and ran deterministic tests and optimization campaigns;
- identified invalid assumptions and withdrew failed routes;
- interpreted certificates, traces and counterexamples;
- proposed narrower obligations and continued from saved checkpoints.

The human operator supplied the goal, constraints, laptop access, review and
continuation decisions. The safe short form is therefore:

> Human-directed, AI-executed research loop with formal checkpoints.

This wording describes execution and provenance. It does not promote finite
evidence into a proof; expert mathematical review remains required.
