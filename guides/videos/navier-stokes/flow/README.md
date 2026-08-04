# Google Flow workflow

Recommended route: **Veo 3.1 Lite at 720p, one output at a time**. Use Fast only
when a Lite take has the correct composition but insufficient motion quality.
Quality is an optional hero-shot upgrade after the complete mobile rough cut,
not a required pass for every shot.

For literal click-by-click instructions, use `RUNBOOK.md`. This
route uses only the standard prompt box and Frames; Agent, Storyboard Studio,
Ingredients and scene automation are not required.

## One-time setup

1. Create a project named `Navier Stokes — Falsify Loop — 9x16`.
2. Use the standard prompt box; turn **Agent** off for the controlled shot pass.
3. Click the model name, choose **Video**, then choose **Veo 3.1 Lite**.
4. Set aspect ratio to **9:16**, resolution to **720p** when exposed, output
   count to **1**, and the nearest supported
   duration listed in `PROMPTS.md`.
5. Drag the shot's keyframe into **Add start frame**.
6. Paste the shared continuity block, one shot block and the negative block into
   the prompt box.
7. Generate once, inspect, and stop when the shot passes. Make one targeted Lite
   retry for a concrete defect. Escalate only that shot to Fast when necessary.

S07 is an editor-held final frame from S06, not a seventh model generation.

Flow feature support changes by model. If the selected Quality model does not
offer a requested frame option, keep the start frame and textual continuity
block; do not switch off 9:16. Ingredients are not used in this workflow.

## Where each text goes

| Flow control | Paste / choose |
| --- | --- |
| Project name | `Navier Stokes — Falsify Loop — 9x16` |
| Mode | Standard prompt box → `Video` |
| Scout model | `Veo 3.1 Lite` |
| Selective upgrade | `Veo 3.1 Fast`; Quality only for an approved hero shot |
| Aspect ratio | `9:16` |
| Resolution | `720p` for all scouts |
| Outputs | 1 per decision cycle |
| Start frame | The reference specified per shot |
| Main prompt | `PROMPTS.md` → continuity + shot + negative |
| S07 overlay copy | `../STORYBOARD.md` → S07; add in the editor only |

## Audio

If Flow offers generated audio, request only subtle glassy water ambience and
no speech. The final English voiceover, score and rhythmic freezes are added
in post so the six shots share one mix.

## Selection rule

Prefer exact composition and filament continuity over spectacular motion. A
slightly quieter shot is easier to speed-ramp than a beautiful take whose cube
warps or whose filament becomes several objects.
