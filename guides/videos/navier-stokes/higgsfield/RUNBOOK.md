# Higgsfield — button-by-button, basic mode

Use this primarily to repair camera problems in S02, S05, or S06. Standard
image-to-video generation is sufficient; camera direction is already written
into the prompts.

## Do not use

- AI Director or automatic storyboards;
- Cast or Soul ID;
- Motion Control or a motion reference;
- a speaking character, lip sync, or automatic narration;
- multi-scene generation.

If Higgsfield places the normal video panel inside a page named `Cinema Studio`,
you may enter that page, but use only **Video / Image to Video**. Director and
character tools are unnecessary.

## Initial setup

1. Open https://higgsfield.ai/ai/video and sign in.
2. Open **Video** or **Image to Video**.
3. Choose `Kling 3.0 Standard`.
4. Set aspect ratio to `9:16`.
5. Turn Audio/Sound off.
6. Set resolution to `720p` and output count to `1`.

## For every shot

1. Upload the PNG listed below into the image/start-frame field.
2. Open `PROMPTS.md`.
3. Paste into the main prompt field:
   - first, the `Shared continuity` fenced block;
   - immediately after it, the selected `S0X` fenced block.
4. Paste `Shared negative` into a separate negative-prompt field when one
   exists. Otherwise append it to the main prompt.
5. Set the listed duration.
6. If the base panel exposes separate lens/camera fields, use the values below.
   If it does not, do not search for a plugin; the prompt contains the same
   information.
7. Click Generate.
8. Download the result into `../renders/higgsfield/`.
9. Review it with `../PRODUCTION_WORKFLOW.md` before generating again.

| Shot | Start image | Length | Camera/lens if separately available |
| --- | --- | ---: | --- |
| S01 | `../assets/keyframe-01-world.png` | 5 s | 50 mm, slow dolly in |
| S02 | `../assets/keyframe-01-world.png` | 5–6 s | 24 mm, push in / tracking |
| S03 | `../assets/keyframe-03-checkpoint.png` | 5 s | 85 mm macro, 30° orbit |
| S04 | `../assets/keyframe-01-world.png` | 8 s | 35 mm, fast fly-through |
| S05 | `../assets/keyframe-03-checkpoint.png` | 8 s | 100 mm macro, locked |
| S06 | `../assets/keyframe-06-frontier.png` | 8 s | 24 mm, slow dolly out |

## Revision rule

If the scout passes, stop generating and place it in the rough cut. If it
fails, change only the failed camera sentence and make one more 720p Standard
attempt. Do not add AI Director, storyboard, or Motion Control. Upgrade quality
only if the complete phone rough cut reveals a visible defect.
