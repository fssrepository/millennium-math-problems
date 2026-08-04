# Navier–Stokes short-video production bundle

Status: production guide for a conceptual 9:16 explainer, not a numerical simulation export  
Generation target: 720 × 1280 first; optional upscale only after picture lock
Timeline target: 1080 × 1920, 30 fps, 29.7 seconds, English narration and captions

## Recommended route

For this exact film, the strongest workflow is:

1. Build the complete 29.7-second still-image animatic before spending credits.
2. Test only S02, S05, and S06 at 720p with one **Veo 3.1 Lite** output each.
3. If those pass, generate S01, S03, and S04 one at a time with Lite.
4. Use **Veo 3.1 Fast** only for a shot that fails visually at Lite quality.
5. Use **Veo 3.1 Quality** only for one or two hero shots if a phone review
   shows a visible need; it is not a mandatory final pass.
6. If a camera still fails after one targeted retry, remake only that shot with
   **Higgsfield or OpenArt basic image-to-video + Kling 3.0 Standard at 720p**.
7. Add every title, caption, arrow, freeze, time-ramp and counter in a normal
   editor such as DaVinci Resolve, Premiere or CapCut. Do not ask a video model
   to draw the English text.

Why this split: the film has no human character, so Flow/Veo's atmosphere and
scene fidelity matter more than face performance. Higgsfield/Kling is the
retake path for camera moves. OpenArt is an optional basic image-to-video
alternative. No Agent, Storyboard Studio, Smart Shot, Director, Motion Control,
character system or third-party plugin is required.

## Model choice at a glance

| Platform | Cheapest useful test | Selective upgrade | Best use here |
| --- | --- | --- | --- |
| Google Flow | Veo 3.1 Lite, 720p, one output | Fast only when needed; Quality only for a visible hero-shot gain | Main fluid footage |
| Higgsfield | Kling 3.0 Standard image-to-video, 720p, one output | Higher mode only after a successful scout | Camera rescue for S02/S05/S06 |
| OpenArt | Kling 3.0 Standard image-to-video, 720p, one output | Higher quality only after a successful scout | Complete fallback or alternate composition |

The recommendation is dated 2026-08-04. Model menus change, so the platform
guides include a fallback rule: preserve the supplied reference frames and
choose the newest model that supports 9:16 plus image/start-frame guidance.

## What is in this bundle

- `00_START_HERE.md` — the first file to open; a plain-language folder map
  and seven-line workflow.
- `PRODUCTION_WORKFLOW.md` — the manual generate → inspect → retry → accept
  workflow with a quality gate after every shot.
- `REVIEW_SHEET.csv` — take-by-take acceptance log.
- `STORYBOARD.md` — the complete 29.7-second shot plan.
- `TIMELINE_EVIDENCE.md` — audited T+ milestones, exact screen mapping and the
  honest next-step estimate.
- `CREATIVE_BIBLE.md` — world, recurring visual “characters,” palette,
  continuity and master prompts.
- `EDIT_PLAN.md` — exact time-remaps, freeze beats, overlays, arrows, sound and
  export settings.
- `COPY_PASTE_MAP.md` — the shortest field-by-field production checklist.
- `EDITOR_OVERLAYS_EN.md` — all English screen copy grouped by timecode.
- `MANUAL_COPY_EN.md` — concise, non-specialist English guide copy.
- `SOCIAL_COPY_EN.md` — title, post copy, pinned clarification and tags.
- `narration/en.md` — final English voiceover and delivery prompt.
- `captions/navier-stokes-en.srt` — timed subtitles.
- `shot-manifest.csv` — production tracking sheet.
- `openart/`, `flow/`, `higgsfield/` — exact UI workflow and copy/paste prompts
  for each platform.
- Each platform contains `RUNBOOK.md` with literal field-by-field steps and
  an explicit list of optional tools not to use.
- `flow/READY_TO_PASTE.md` — fully assembled one-block prompts for the default
  Flow workflow.
- `renders/` — prepared local drop zones for drafts, finals, audio and edits.
- `assets/` — three consistent 9:16 keyframes and their provenance.
- `sources/RESEARCH_NOTES.md` — repository/session evidence and current product
  sources.
- `DELIVERY_QC.md` — final fact, visual, audio and platform checks.

## Fastest way to make it

1. Open `00_START_HERE.md`.
2. Follow `PRODUCTION_WORKFLOW.md` and `flow/RUNBOOK.md`.
3. Generate and approve shots in scout order S02 → S05 → S06 → S01 → S03 →
   S04, one output at a time.
4. Assemble them according to `STORYBOARD.md` and `EDIT_PLAN.md`.
5. Import the SRT, record or synthesize `narration/en.md`, then run
   `DELIVERY_QC.md`.

## Non-negotiable scientific wording

- The Millennium problem is **not solved**.
- The footage is a **concept visualization**, not the output of the C++ run.
- `F000`–`F015` means **16 rejected lemma candidates**.
- A finite search can falsify a candidate or let it survive the tested set; it
  cannot prove the cutoff-independent L4 statement.
- The demonstrated result is a fast, reproducible lemma-falsification loop,
  one proved partial far-tail estimate, and a still-open L4 closure.
- The `~1–3 days` estimate applies only to the next PNT-12 decision, not to a
  proof and not to the Clay problem.
- Protas/Tao is a methodological comparison, not provenance of a combined
  proof; use the exact wording in `STORYBOARD.md`.

## The short version of the story

`Why → fluid problem → method family → fast early failures → slower deep
optimization → surviving open frontier → honest estimate.`

From 00:12.40 to 00:25.20 the on-screen research clock is a linear compression
of 33 h 28 min. Early failures arrive quickly; later frames breathe longer as
the searches move from scaling and triads to exact gradients and K12.
