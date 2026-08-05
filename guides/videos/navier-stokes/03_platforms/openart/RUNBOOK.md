# OpenArt — button-by-button, basic mode

This is a complete alternative workflow using the standard **Image to Video**
surface.

## Do not use

- One Click Story;
- Smart Shot;
- OpenArt Director;
- Motion Control or a motion-reference video;
- AI Character, lip sync, or a built-in narrator.

## Initial setup

1. Open https://openart.ai/ai-video-generator/ and sign in.
2. Open **Create Video**.
3. Select **Image to Video**. A field labeled only `Start Frame` is also valid.
4. Choose `Kling 3.0 Standard` image-to-video. Do not choose `Motion
   Control`; this film has no human motion reference.
5. Set aspect ratio to `9:16`.
6. Set resolution to `720p` and output count to `1`.
7. Turn audio/voice off or choose ambient only.

## For every shot

1. Upload the start image listed below.
2. Open `PROMPTS.md`.
3. Paste into the main **Prompt** field in this order:
   - the complete `Shared continuity` fenced block;
   - the selected `S0X` fenced block.
4. If a **Negative prompt** field exists, paste `Shared negative` there. If it
   does not, append that block to the main prompt.
5. Set the shot duration.
6. Click **Generate**.
7. Download the result into `../../renders/openart/`.
8. Review it with `../../02_production/PRODUCTION_WORKFLOW.md` before generating again.

| Shot | Start image | Length |
| --- | --- | ---: |
| S01 | `../../assets/keyframe-01-world.png` | 5 s |
| S02 | `../../assets/keyframe-01-world.png` | 6 s or nearest available |
| S03 | `../../assets/keyframe-03-checkpoint.png` | 5 s |
| S04 | `../../assets/keyframe-01-world.png` | 8–10 s |
| S05 | `../../assets/keyframe-03-checkpoint.png` | 8–10 s |
| S06 | `../../assets/keyframe-06-frontier.png` | 8–10 s |

## After reviewing the scout

1. If the scout passes, stop generating this shot and use it in the rough cut.
2. If it fails, record the defect and seed in `../../02_production/REVIEW_SHEET.csv`.
3. Preserve the image and seed, change one prompt sentence, and make one more
   720p Standard attempt.
4. Upgrade quality only after the full phone rough cut demonstrates a visible
   need in this exact shot.

If the camera remains unstable after one targeted retry, do not open Smart
Shot. Move only that shot back to Flow Fast or to Higgsfield Standard—choose
one, not both at once.
