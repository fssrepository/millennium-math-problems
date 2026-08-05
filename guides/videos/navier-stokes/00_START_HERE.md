# OPEN THIS FILE FIRST

The default workflow uses Google Flow. It requires no plugin, Agent, or
Storyboard Studio. Free Flow access is account- and region-dependent; if the
site does not let the account reach the generation screen, do not subscribe
automatically. Use the comparison in `02_production/PLATFORM_COST_GUIDE.md` and select one
paid fallback instead.

Open: https://flow.google

## Where everything is

| What you need | File or directory |
| --- | --- |
| Start-frame images to upload | `assets/` |
| Story, evidence and narration | `01_story/` |
| Platform, model and spending decision | `02_production/PLATFORM_COST_GUIDE.md` |
| Current Flow controls and exact prompts | `03_platforms/flow/CHAIN_32S_RUNBOOK.md` |
| Legacy independent-shot workflow | `03_platforms/flow/RUNBOOK.md` |
| Checks to run after every generation | `02_production/PRODUCTION_WORKFLOW.md` |
| Take-selection log | `02_production/REVIEW_SHEET.csv` |
| Current production state and credit ledger | `02_production/PRODUCTION_LOG.md` |
| Current ten-frame visual plan | `assets/previs-10/README.md` |
| Evidence-driven PNT-13 failure shot | `01_story/FALSIFICATION_VISUAL.md` |
| Reproducible local video repairs | `scripts/` |
| Rejected/superseded media decisions | `02_production/REVIEW_SHEET.csv` and `02_production/PRODUCTION_LOG.md` |
| Local animatic build | `scripts/build_animatic.sh` (review export is not a backup to commit) |
| Exact clip order and duration | `02_production/shot-manifest.csv` |
| Speed ramps and freeze frames | `02_production/EDIT_PLAN.md` |
| English on-screen copy | `04_editorial/EDITOR_OVERLAYS_EN.md` |
| Subtitle file | `04_editorial/captions/navier-stokes-en.srt` |
| Narration and TTS direction | `01_story/narration/en.md` |
| Downloaded clip destinations | `renders/` |

Expanded Agent fallback: after the Agent describes the prepared job, type
`Approve` as a separate prompt message. Then inspect and click `Approve` on the
built-in cost card. The typed word does not spend credits; the clicked option
does. Never select `Approve, do not ask again`.

## Current review route

Read `01_story/FALSIFICATION_VISUAL.md`, then inspect
`renders/edit/navier-stokes-10-frame-story-previs-v02-30s.mp4` at phone size.
Do not approve V01 until its continuous camera move, tail inspection and final
application morph are accepted. The paid route then starts fresh: one
eight-second Lite video followed by three Extend actions on the same video.

Run `scripts/build_animatic.sh`, then watch its local 29.70-second output at
phone size. It labels every shot as either a provisional Flow plate or a still
placeholder and uses raw S02/S05 motion without synthetic zoom, slowdown or
freeze. Do not keep superseded review exports in Git; the review decision goes
in `02_production/REVIEW_SHEET.csv`.

## The Flow route after animatic review

1. Read `02_production/PLATFORM_COST_GUIDE.md`, then open the Google Flow website. If the
   generation screen is available, continue with
   `03_platforms/flow/CHAIN_32S_RUNBOOK.md`. Otherwise
   stop this list, choose either Higgsfield or OpenArt from the cost table, and
   follow that platform's `RUNBOOK.md` from top to bottom.
2. Generate V01 once with Veo 3.1 Lite: 9:16, 720p, eight seconds, one output.
3. Inspect V01, then use Extend on that same accepted video for V02, V03 and
   V04. Never create a fresh scene at a slot boundary.
4. Require a cost of no more than 10 credits for each action. Four actions use
   40 credits; keep the last 10 as one targeted retry reserve.
5. Generate the locked instrumental score in V01 and continue the same cue
   through V02–V04; do not add a separate music track. Reject an audible
   restart, click or style change at any Extend boundary.
6. Add all guides, labels, status graphics, captions and narration in post.

If a camera move still fails after one targeted revision, move only that shot
to `03_platforms/higgsfield/RUNBOOK.md`. The complete OpenArt alternative is documented in
`03_platforms/openart/RUNBOOK.md`.

Alternative interfaces:

- OpenArt video: https://openart.ai/ai-video-generator/
- Higgsfield video: https://higgsfield.ai/ai/video
