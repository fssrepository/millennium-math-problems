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
| Superseded continuity-only morph previs | `renders/edit/navier-stokes-10-frame-morph-previs-v01-30s.mp4` |
| Reproducible local video repairs | `scripts/` |
| Current 29.70-second animatic | `renders/edit/navier-stokes-animatic-v02-raw-plates-720x1280.mp4` |
| Current untrimmed Flow review | `renders/edit/navier-stokes-raw-flow-review-v01-44p2s.mp4` |
| Current four-source Flow review | `renders/edit/navier-stokes-four-flow-clips-v01-32s.mp4` |
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

Read `FALSIFICATION_VISUAL.md`, then use
`renders/edit/navier-stokes-10-frame-morph-previs-v01-30s.mp4` only to review
the broad camera continuity. The v01 morph predates the explicit expected
decay → surviving tail → saved PNT-13 rejection → PNT-12 open sequence and is
not story-approved. Do not paste or approve C01/C02/C03 until a v02 preview
contains those four visible state changes. The first Flow video remains;
regeneration restarts from the second video piece.

Run `scripts/build_animatic.sh`, then watch
`renders/edit/navier-stokes-animatic-v01-720x1280.mp4`. It is exactly 29.70
seconds and labels every shot as either a provisional Flow plate or a still
placeholder. The current build uses raw S02 and S05 motion without synthetic
zoom, slowdown or freeze. Do not generate another shot until this complete
timeline has been reviewed at phone size.

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
