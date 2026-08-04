# OPEN THIS FILE FIRST

The default workflow uses Google Flow. It requires no plugin, Agent, or
Storyboard Studio.

Open: https://flow.google

## Where everything is

| What you need | File or directory |
| --- | --- |
| Start-frame images to upload | `assets/` |
| Exact Flow button sequence | `flow/RUNBOOK.md` |
| Complete Flow prompts ready to paste | `flow/READY_TO_PASTE.md` |
| Checks to run after every generation | `PRODUCTION_WORKFLOW.md` |
| Take-selection log | `REVIEW_SHEET.csv` |
| Exact clip order and duration | `shot-manifest.csv` |
| Speed ramps and freeze frames | `EDIT_PLAN.md` |
| English on-screen copy | `EDITOR_OVERLAYS_EN.md` |
| Subtitle file | `captions/navier-stokes-en.srt` |
| Narration and TTS direction | `narration/en.md` |
| Downloaded clip destinations | `renders/` |

## The whole workflow in seven sentences

1. Open `flow/RUNBOOK.md` and the Google Flow website.
2. In Flow select: standard prompt box → Video → Frames → 9:16.
3. Upload the image named for that shot from `assets/`.
4. Paste the complete matching block from `flow/READY_TO_PASTE.md`.
5. Generate four drafts, download them, and inspect them with
   `PRODUCTION_WORKFLOW.md`.
6. Make two Quality versions of the accepted draft, then repeat for S01–S06.
7. Assemble the six accepted clips with `EDIT_PLAN.md`, then add the copy from
   `EDITOR_OVERLAYS_EN.md`, the SRT file, and the narration.

If a camera move still fails after two targeted revisions, move only that shot
to `higgsfield/RUNBOOK.md`. The complete OpenArt alternative is documented in
`openart/RUNBOOK.md`.

Alternative interfaces:

- OpenArt video: https://openart.ai/ai-video-generator/
- Higgsfield video: https://higgsfield.ai/ai/video
