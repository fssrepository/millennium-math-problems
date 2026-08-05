# Google Flow workflow

> **CURRENT ROUTE:** use `CHAIN_32S_RUNBOOK.md`. It creates one continuous
> video as V01 plus three Veo 3.1 Lite Extends: four 8-second slots, 40 credits,
> and a 10-credit retry reserve. The independent-shot files below remain only
> as production history.

Recommended route: **Veo 3.1 Lite at 720p, one output at a time**. Use Fast only
when a Lite take has the correct composition but insufficient motion quality.
Quality is an optional hero-shot upgrade after the complete mobile rough cut,
not a required pass for every shot.

For current click-by-click instructions and paste-ready prompts, use
`CHAIN_32S_RUNBOOK.md`.

If Flow opens the expanded Agent fallback used in the supervised pass, follow
`AGENT_READY_TO_PASTE.md`. After the Agent describes the prepared job, send the
single word `Approve`; then inspect and click `Approve` on Flow's built-in cost
card. Typing the word does not spend credits, while clicking the built-in
option does. Never choose `Approve, do not ask again`.

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
| S07 overlay copy | `../../01_story/STORYBOARD.md` → S07; add in the editor only |

## Audio

For the current V01–V04 Extend route, generate the musical bed in V01 and carry
the exact same cue through every Extend. The locked style is restrained premium
science-documentary electronica with a mid-tempo low pulse, glassy-water
texture, sparse dry ticks and a soft airy whoosh. No speech, singing or lyrics.
The full boundary-safe wording and audio rejection checks are in
`CHAIN_32S_RUNBOOK.md`.

This removes the need to source or add a separate music track. English
voiceover, captions and factual graphics remain deterministic editorial
elements; do not ask Veo to generate them.

## Selection rule

Prefer exact composition and filament continuity over spectacular motion. A
slightly quieter shot is easier to speed-ramp than a beautiful take whose cube
warps or whose filament becomes several objects.
