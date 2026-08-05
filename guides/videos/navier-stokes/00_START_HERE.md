# OPEN THIS FILE FIRST

The default workflow uses Google Flow. It requires no plugin, Agent, or
Storyboard Studio. Free Flow access is account- and region-dependent; if the
site does not let the account reach the generation screen, do not subscribe
automatically. Use the comparison in `PLATFORM_COST_GUIDE.md` and select one
paid fallback instead.

Open: https://flow.google

## Where everything is

| What you need | File or directory |
| --- | --- |
| Start-frame images to upload | `assets/` |
| Platform, model and spending decision | `PLATFORM_COST_GUIDE.md` |
| Exact Flow button sequence | `flow/RUNBOOK.md` |
| Complete Flow prompts ready to paste | `flow/READY_TO_PASTE.md` |
| Current exact prompt for the expanded Agent UI | `flow/AGENT_READY_TO_PASTE.md` |
| Checks to run after every generation | `PRODUCTION_WORKFLOW.md` |
| Take-selection log | `REVIEW_SHEET.csv` |
| Current production state and credit ledger | `PRODUCTION_LOG.md` |
| Current ten-frame visual plan | `assets/previs-10/README.md` |
| Evidence-driven PNT-13 failure shot | `FALSIFICATION_VISUAL.md` |
| Reproducible local video repairs | `scripts/` |
| Rejected/superseded media decisions | `REVIEW_SHEET.csv` and `PRODUCTION_LOG.md` |
| Local animatic build | `scripts/build_animatic.sh` (review export is not a backup to commit) |
| Exact clip order and duration | `shot-manifest.csv` |
| Speed ramps and freeze frames | `EDIT_PLAN.md` |
| English on-screen copy | `EDITOR_OVERLAYS_EN.md` |
| Subtitle file | `captions/navier-stokes-en.srt` |
| Narration and TTS direction | `narration/en.md` |
| Downloaded clip destinations | `renders/` |

Expanded Agent fallback: after the Agent describes the prepared job, type
`Approve` as a separate prompt message. Then inspect and click `Approve` on the
built-in cost card. The typed word does not spend credits; the clicked option
does. Never select `Approve, do not ask again`.

## Current review route

Read `FALSIFICATION_VISUAL.md`, then inspect the ten-frame plan under
`assets/previs-10/`. The superseded morph-previs binary was removed during media
cleanup; its decision remains in `PRODUCTION_LOG.md`. Do not paste or approve
C01/C02/C03 until a locally built v02 preview contains expected decay →
surviving tail → saved PNT-13 rejection → PNT-12 open. The first Flow video
remains; regeneration restarts from the second video piece.

Run `scripts/build_animatic.sh`, then watch its local 29.70-second output at
phone size. It labels every shot as either a provisional Flow plate or a still
placeholder and uses raw S02/S05 motion without synthetic zoom, slowdown or
freeze. Do not keep superseded review exports in Git; the review decision goes
in `REVIEW_SHEET.csv`.

## The Flow route after animatic review

1. Read `PLATFORM_COST_GUIDE.md`, then open the Google Flow website. If the
   generation screen is available, continue with `flow/RUNBOOK.md`. Otherwise
   stop this list, choose either Higgsfield or OpenArt from the cost table, and
   follow that platform's `RUNBOOK.md` from top to bottom.
2. In Flow select: standard prompt box → Video → Frames → 9:16.
3. Upload the image named for that shot from `assets/`.
4. Paste the complete matching block from `flow/READY_TO_PASTE.md`.
5. Generate one selected 720p Lite scout, download it, and inspect it with
   `PRODUCTION_WORKFLOW.md`.
6. Continue only if it passes; otherwise stop that shot after one targeted
   retry. Use the animatic review—not the old scout order—to choose the next
   highest-impact missing plate.
7. Replace placeholders in the animatic one at a time. Upgrade only visibly
   weak shots, then retain the copy from
   `EDITOR_OVERLAYS_EN.md`, the SRT file, and the narration.

If a camera move still fails after one targeted revision, move only that shot
to `higgsfield/RUNBOOK.md`. The complete OpenArt alternative is documented in
`openart/RUNBOOK.md`.

Alternative interfaces:

- OpenArt video: https://openart.ai/ai-video-generator/
- Higgsfield video: https://higgsfield.ai/ai/video
