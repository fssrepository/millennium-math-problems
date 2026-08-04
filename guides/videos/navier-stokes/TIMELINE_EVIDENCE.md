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
| 2026-08-02 17:19:22 | T+33:28 | 00:25.20 | PNT-12 finite stress campaign reported | `SURVIVED THE TEST` |

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

## Provenance

Primary local record:

`/home/raxim/.codex/archived_sessions/rollout-2026-08-01T09-46-19-019fbc49-c7f3-7051-8627-e920169afd1b.jsonl`

Rejected-lemma ledger:

`../../../navier-stokes/proof/failed_lemmas.tsv`

The private session archive is referenced for auditability but is not copied
into this production bundle.
