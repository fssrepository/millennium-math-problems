# Higgsfield workflow

Recommended route: **Cinema Studio 3.5 + Kling 3.0**. Kling is the primary
model because this film needs macro stability and deliberate camera movement,
not a speaking human. Use Veo 3.1 inside Higgsfield only as a final-fidelity
alternate for S01 or S06 if the chosen plan exposes it with the same reference
controls.

## One-time setup

1. Open **Cinema Studio** and choose **9:16**.
2. Choose **Kling 3.0** as the generation model.
3. Upload the shot's keyframe as image/start-frame reference.
4. Set the director controls from the table below. If an exact preset name is
   unavailable, choose its nearest neutral equivalent and keep the written
   prompt unchanged.
5. Paste `PROMPTS.md` → shared continuity + one shot prompt into the scene
   prompt. Put the shared negative block in the negative field if present;
   otherwise append it.
6. Render 3 variants. Change only one control between revisions.

## Director controls

| Shot | Lens | Camera move | Focus | Lighting | Duration |
| --- | --- | --- | --- | --- | --- |
| S01 | 50 mm | Slow dolly in | Deep enough for full cube | Low-key cyan rim | 5 s |
| S02 | 24 mm | Fast push in / tracking | Red bend at end | Cyan with restrained coral | 5 s |
| S03 | 85 mm macro | 30° orbit | Red bend | Macro cyan/coral | 5 s |
| S04 | 35 mm | Fast forward fly-through | Shell crossings | Dark science-documentary | 8 s |
| S05 | 100 mm macro | Locked / micro-slider | Scan intersection | Amber scan accent | 8 s |
| S06 | 24 mm | Slow precise dolly out | Cube then shell stack | Open cyan gap | 8 s |

Suggested general settings where available: neutral cinematic genre, high
prompt adherence, low stylization, realistic motion, restrained grain, no
automatic dialogue.

## Camera-effect rule

Higgsfield supplies the clean camera plate. The precise speed ramps, hard
freezes and repeated S05 stop/restart rhythm still belong in the editor. This
keeps the explanation frame-accurate and prevents the model from inventing
text or changing the scene at every “restart.”

## Selection rule

Pick the take whose endpoint composition matches the next keyframe. Camera
spectacle is secondary to a stable red filament and unchanged glass cube.
S07 is the final held frame of S06 with editor typography; do not generate it.
