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
5. Generate one 720p Lite scout, download it, and inspect it with
   `PRODUCTION_WORKFLOW.md`.
6. Continue only if it passes; otherwise make one targeted retry. Process shots
   in cost-risk order S02 → S05 → S06 → S01 → S03 → S04.
7. Assemble the six accepted scouts with `EDIT_PLAN.md`. Upgrade only visibly
   weak shots, then add the copy from
   `EDITOR_OVERLAYS_EN.md`, the SRT file, and the narration.

If a camera move still fails after one targeted revision, move only that shot
to `higgsfield/RUNBOOK.md`. The complete OpenArt alternative is documented in
`openart/RUNBOOK.md`.

Alternative interfaces:

- OpenArt video: https://openart.ai/ai-video-generator/
- Higgsfield video: https://higgsfield.ai/ai/video
