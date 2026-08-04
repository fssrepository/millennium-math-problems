# Google Flow — button-by-button

This is the primary route. It uses only Flow's standard prompt box and built-in
**Frames** image-to-video mode.

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
6. Click the model name beside the prompt box. It may default to an image model
   such as Nano Banana Pro.
7. In the selector, click **Video**.
8. For video input type, choose **Frames**.
9. Drag the shot's PNG onto `+ Add start frame`.
10. Leave `+ Add end frame` empty.
11. Select:
    - aspect ratio: `9:16`;
    - model: `Veo 3.1 Fast`;
    - outputs: `4`;
    - generation length: the nearest available value in the table below.
12. Choose ambient/no speech for audio, or turn sound off when a separate
    switch is available.

If a label changes slightly, the path is still: **standard prompt → Video →
Frames → Add start frame**.

## Process shots in this order

| Order | Shot | Start frame | Length | Prompt to paste |
| ---: | --- | --- | ---: | --- |
| 1 | S01 | `../assets/keyframe-01-world.png` | 5 s | `READY_TO_PASTE.md` → complete S01 block |
| 2 | S02 | `../assets/keyframe-01-world.png` | 6 s | `READY_TO_PASTE.md` → complete S02 block |
| 3 | S03 | `../assets/keyframe-03-checkpoint.png` | 5 s | `READY_TO_PASTE.md` → complete S03 block |
| 4 | S04 | `../assets/keyframe-01-world.png` | 8 s | `READY_TO_PASTE.md` → complete S04 block |
| 5 | S05 | `../assets/keyframe-03-checkpoint.png` | 8 s | `READY_TO_PASTE.md` → complete S05 block |
| 6 | S06 | `../assets/keyframe-06-frontier.png` | 8 s | `READY_TO_PASTE.md` → complete S06 block |

Example for S01:

1. Drag `keyframe-01-world.png` onto `+ Add start frame`.
2. Open `READY_TO_PASTE.md`.
3. Find `S01 — 5 s`.
4. Select the entire English block between the three backticks.
5. Paste it into Flow's main prompt box.
6. Confirm: Video, Frames, 9:16, Veo 3.1 Fast, 4 outputs, 5 s.
7. Click **Generate**. Some localizations may display `Generate Image`; the
   selected Video mode is what matters.
8. Wait for all four outputs.
9. Download them into `../renders/flow/drafts/` as
   `S01_flow_fast_take01.mp4` through `take04.mp4`.
10. Run the review gate in `../PRODUCTION_WORKFLOW.md`.

Repeat with the next row. For every new shot, remove the previous start frame
and prompt before loading the next pair.

## After accepting a Fast draft

1. Do not change the start frame or prompt.
2. Click the model name again.
3. Switch to `Veo 3.1 Quality`.
4. Set outputs to `2`.
5. Keep the same duration and 9:16 ratio.
6. Click Generate.
7. Download both results into `../renders/flow/finals/`.
8. Watch both at normal speed and 0.5× speed.
9. Append `_SELECTED.mp4` to the winner's filename.

## Do not proceed when

- the cube deforms;
- more than one coral filament appears;
- unreadable AI writing or an equation appears;
- the camera misses the requested endpoint;
- fewer than 8–12 clean frames exist at either end.

Do not open Storyboard Studio or Agent to fix these problems. Change one
sentence using the troubleshooting table in `../PRODUCTION_WORKFLOW.md`, then
generate again.
