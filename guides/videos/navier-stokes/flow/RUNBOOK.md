# Google Flow — cost-efficient button-by-button workflow

This is the primary route. It uses only Flow's standard prompt box and built-in
**Frames** image-to-video mode. Generate one 720p result, inspect it, and spend
again only when the result identifies a concrete next action.

## Do not open

- Agent;
- Storyboard Studio or a scene builder;
- Ingredients;
- character, avatar, or voice references;
- Extend or automatic story generation.

## Initial setup

1. Open https://flow.google and sign in.
2. Create a new project.
3. Name it `Navier Stokes — Lemma Stress Test — 9x16`.
4. Locate the large prompt box near the center or bottom of the page.
5. If **Agent** is enabled, turn it off.
6. Click the model name beside the prompt box.
7. Select **Video**.
8. Select **Frames** as the input type.
9. Drag the shot's PNG onto `+ Add start frame`.
10. Leave `+ Add end frame` empty.
11. Select:
    - aspect ratio: `9:16`;
    - model: `Veo 3.1 Lite`;
    - resolution: `720p` when a resolution control is shown;
    - outputs: `1`;
    - generation length: the nearest available value listed below.
12. Choose ambient/no speech for audio, or turn sound off.

If a label changes slightly, the path is still: **standard prompt → Video →
Frames → Add start frame**. If Lite does not expose start-frame input in your
account, use Fast with the same one-output rule.

## Scout in this order

The order minimizes wasted generations by testing the hardest camera behavior
first. It is not the final edit order.

| Scout order | Shot | Start frame | Length | Prompt to paste |
| ---: | --- | --- | ---: | --- |
| 1 | S02 | `../assets/keyframe-01-world.png` | 6 s | `READY_TO_PASTE.md` → complete S02 block |
| 2 | S05 | `../assets/keyframe-03-checkpoint.png` | 8 s | `READY_TO_PASTE.md` → complete S05 block |
| 3 | S06 | `../assets/keyframe-06-frontier.png` | 8 s | `READY_TO_PASTE.md` → complete S06 block |
| 4 | S01 | `../assets/keyframe-01-world.png` | 5 s | `READY_TO_PASTE.md` → complete S01 block |
| 5 | S03 | `../assets/keyframe-03-checkpoint.png` | 5 s | `READY_TO_PASTE.md` → complete S03 block |
| 6 | S04 | `../assets/keyframe-01-world.png` | 8 s | `READY_TO_PASTE.md` → complete S04 block |

## Exact process for one shot

Example: S02.

1. Drag `keyframe-01-world.png` onto `+ Add start frame`.
2. Open `READY_TO_PASTE.md` and find `S02 — 6 s`.
3. Copy the complete block between the three backticks.
4. Paste it into Flow's main prompt box.
5. Confirm: Video, Frames, 9:16, Lite, 720p, 1 output, 6 s.
6. Read the displayed credit charge before clicking. It must describe one
   generation. If the interface forces two outputs, cancel and recheck the
   output setting instead of approving an accidental double generation.
7. Click **Generate**.
8. Download the result into `../renders/flow/drafts/` as
   `S02_flow_lite_720p_take01.mp4`.
9. Watch it at normal speed and 0.5× speed.
10. Complete the S02 row in `../REVIEW_SHEET.csv`.
11. If it passes, stop generating S02 and move to S05.
12. If it fails, choose one repair from `../PRODUCTION_WORKFLOW.md`, change only
    that sentence, and generate one more Lite output.
13. If the same defect remains, stop using Lite for this shot. Try one Fast
    output or one Kling Standard fallback.

Repeat with the next row. Remove the previous start frame and prompt before
loading the next shot.

## When to use Fast

Use `Veo 3.1 Fast` only when:

- the Lite composition and endpoint are correct;
- a motion or detail defect remains visible at phone size;
- the full 720p rough cut shows that this particular shot needs improvement.

Keep 720p, one output, the same start frame, and the accepted prompt. Append
`_fast_take01.mp4` to the filename.

## When to use Quality

Quality is optional. Consider it only after the complete rough cut passes, and
only for S02 or S06 when an A/B test shows a visible phone-sized improvement.
Generate one output, trim it in the editor, and do not run Quality across all
six shots by default.

## Do not proceed when

- the cube deforms;
- more than one coral filament appears;
- unreadable AI writing or an equation appears;
- the camera misses the requested endpoint;
- fewer than 8–12 clean frames exist at either end.

Do not open Storyboard Studio or Agent to fix these problems.
