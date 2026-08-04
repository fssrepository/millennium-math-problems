# Higgsfield workflow

Recommended route: the platform's basic **Image to Video + Kling 3.0 Standard**
surface at 720p, one output at a time, and only for a specific failed Flow shot.
The camera and lens instructions are already in the text prompts. AI Director,
storyboard, Cast, Soul ID and Motion Control are not required. If the normal
video panel appears inside a page named Cinema Studio, use only that base panel.

For literal click-by-click instructions, use `RUNBOOK.md`.

## One-time setup

1. Open **Video / Image to Video** and choose **9:16**.
2. Choose **Kling 3.0** as the generation model.
3. Upload the shot's keyframe as image/start-frame reference.
4. If the base panel exposes lens/camera controls, use the table below. If not,
   do nothing: the same information is already in the prompt.
5. Paste `PROMPTS.md` → shared continuity + one shot prompt into the scene
   prompt. Put the shared negative block in the negative field if present;
   otherwise append it.
6. Render one 720p scout. Inspect it before authorizing another generation.

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
