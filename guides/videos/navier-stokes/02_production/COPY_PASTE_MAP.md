# Copy/paste map

Use this page when producing; it is the shortest route through the bundle.

> **CURRENT HOLD:** do not paste the legacy S05/S06 prompts. First read
> `../01_story/FALSIFICATION_VISUAL.md` and approve a v02 preview containing the PNT-13
> failure and PNT-12 state transition. The platform prompt files are marked
> `HOLD` until that rewrite is complete.

## 1. Pick one video platform

Recommended default: **Google Flow when the account exposes free generation**.
If it does not, compare Higgsfield and OpenArt in `PLATFORM_COST_GUIDE.md` and
pay for only one of them.

For the fewest copy operations, follow `../00_START_HERE.md` and use the complete
blocks in `../03_platforms/flow/READY_TO_PASTE.md`.
For literal button-by-button instructions, open `../03_platforms/flow/RUNBOOK.md`.

- Project/settings fields: copy the values from `../03_platforms/flow/README.md` →
  **Where each text goes**.
- Main prompt field: copy `../03_platforms/flow/PROMPTS.md` → **Shared continuity**, then the
  chosen **S01–S06** block, then **Shared negative**.
- Start-frame field: upload the PNG named in that shot heading.
- Generated-audio field, when present: ambient audio only; no speech.

Alternative paths:

- OpenArt: follow `../03_platforms/openart/RUNBOOK.md`, then use
  `../03_platforms/openart/PROMPTS.md` in the same continuity → shot → negative order.
- Higgsfield: follow `../03_platforms/higgsfield/RUNBOOK.md`; use its basic Image to Video
  panel, put continuity + shot into the main prompt and the negative block into
  the negative-prompt field.

Generate only S01–S06. S07 is a held final S06 frame.

## 2. Put these texts in the editor—not in the video model

- Shot titles, event chips and final card: `../01_story/STORYBOARD.md` → **Post overlay**.
- Expected decay, persistent far-shell response and PNT-12 handoff:
  `../01_story/FALSIFICATION_VISUAL.md`.
- Exact T+ placement: `../01_story/TIMELINE_EVIDENCE.md` → **Audited milestones**.
- Speed ramps, freezes, pointer and crash/resume dip: `EDIT_PLAN.md`.
- Subtitles: import `../04_editorial/captions/navier-stokes-en.srt`.
- Voice recording/TTS script: `../01_story/narration/en.md` → **Final script**.
- TTS voice-style field: `../01_story/narration/en.md` → fenced text under
  **Voice / TTS direction**.

Never paste English captions into Flow, OpenArt or Higgsfield. Generate clean
visual plates and add typography in CapCut, DaVinci Resolve or Premiere.

## 3. Model order

1. Flow / Veo 3.1 Lite: one 720p scout at a time.
2. Inspect before spending on the next generation; make at most one targeted
   Lite retry for a specific defect.
3. Use Flow / Veo 3.1 Fast only when a selected Lite take is visibly too weak.
4. Use Quality only for one or two hero shots after the 720p rough cut passes a
   phone review.
5. If S02, S05 or S06 camera control still fails, remake only that shot with
   Higgsfield or OpenArt basic Image to Video + Kling 3.0 Standard at 720p.

This route does not require Agent, Storyboard Studio, Smart Shot, Director,
Motion Control, Cast, Soul ID or any third-party plugin.

The recommendation is a 2026-08-04 snapshot; if menus change, preserve 9:16,
the supplied start frame and the continuity block, then select the newest model
that supports those controls.
