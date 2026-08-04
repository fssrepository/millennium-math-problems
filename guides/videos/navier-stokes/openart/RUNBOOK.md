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
4. Choose the normal `Kling 3.0` image-to-video model. Do not choose `Motion
   Control`; this film has no human motion reference.
5. Set aspect ratio to `9:16`.
6. Set output count to `4` for drafts.
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
7. Download all four results into `../renders/openart/`.
8. Review them with `../PRODUCTION_WORKFLOW.md`.

| Shot | Start image | Length |
| --- | --- | ---: |
| S01 | `../assets/keyframe-01-world.png` | 5 s |
| S02 | `../assets/keyframe-01-world.png` | 6 s or nearest available |
| S03 | `../assets/keyframe-03-checkpoint.png` | 5 s |
| S04 | `../assets/keyframe-01-world.png` | 8–10 s |
| S05 | `../assets/keyframe-03-checkpoint.png` | 8–10 s |
| S06 | `../assets/keyframe-06-frontier.png` | 8–10 s |

## Finalize the accepted draft

1. Record the seed in `../REVIEW_SHEET.csv` when OpenArt exposes one.
2. Preserve the image, prompt, and seed.
3. Set output count to 2 and select higher quality when a quality switch exists.
4. Generate again.
5. Append `_SELECTED.mp4` to the winning filename.

If the camera remains unstable after two targeted retries, do not open Smart
Shot. Move only that shot to Higgsfield's basic image-to-video workflow.
