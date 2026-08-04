# OpenArt workflow

Recommended OpenArt route: use **Kling 3.0 Standard at 720p, one output at a
time**, and only after a Flow scout identifies a reason to switch. Use the basic
9:16 image-to-video/start-frame workflow. Smart Shot,
One Click Story, Director and Motion Control are not required.

For literal click-by-click instructions, use `RUNBOOK.md`.

## One-time setup

1. Open **Create Video**.
2. Select **9:16**.
3. Select **Start + End Frame** or **Image to Video** when available for the
   chosen model.
4. Select **Kling 3.0**. If its exact listing says `Omni` or `Motion Control`,
   use the normal start-frame version for this film; no body-motion reference
   is needed.
5. Upload `../assets/keyframe-01-world.png` as the visual reference for S01,
   S02 and S04; `keyframe-03-checkpoint.png` for S03/S05; and
   `keyframe-06-frontier.png` for S06. S07 uses the final S06 frame as a still
   in the editor; it does not need a separate generation.
6. Paste the shared continuity block, the selected shot prompt, then the shared
   negative block from `PROMPTS.md` into the prompt field.
7. Generate one 720p scout and inspect it. Make one targeted retry only when a
   specific defect is identified.

## Where each text goes

| OpenArt field | Paste / choose |
| --- | --- |
| Aspect ratio | `9:16` |
| Video model | `Kling 3.0 Standard` |
| Start image | The shot reference listed above |
| Prompt | `PROMPTS.md` → shared continuity + one shot block |
| Negative prompt, if shown | `PROMPTS.md` → shared negative block |
| Duration | S01/S03 5 s; S02 6 s; S04–S06 8–10 s |
| Resolution | `720p` for scouts |
| Output count | 1 per decision cycle |

Do not use One Click Story, Smart Shot or Director. This film depends on
editor-accurate freezes, counters and English typography. Do not ask OpenArt
to render the captions.

## Selection rule

Reject a take if the cube deforms, the red filament duplicates, the fluid looks
like smoke, or the camera loses the red filament before the requested end
frame. Keep the same seed if OpenArt exposes it, but the uploaded keyframes are
the primary continuity anchor.
