# Production log

This file persists the decisions made during the supervised production pass.
It is the durable record; chat instructions are not the source of truth.

Last updated: 2026-08-05  
Current next action: rebuild the ten-frame story as v02 with the explicit
PNT-13 falsification sequence; spend no additional Flow credits and do not
finalize prompts before that review

Media cleanup on 2026-08-05 removed rejected, failed-repair, review-only and
superseded MP4 binaries from the working tree. The historical filenames below
remain as an audit trail only; `REVIEW_SHEET.csv` carries the decisions. Git now
retains only active selected sources, accepted derivatives and deliberate final
deliverables. No repository history was rewritten.

Superseded continuity-only planning review:
`../renders/edit/navier-stokes-10-frame-morph-previs-v01-30s.mp4` — 720 × 1280,
24 fps, exactly 30.000 seconds. It cross-morphs ten keyframes at three-second
story intervals. K01–K03 come directly from the accepted first Flow video;
K04–K10 plan the replacement route from the second video piece onward. It does
not yet show expected PNT-13 decay, the measured persistent far-shell response,
the saved rejection or the PNT-12 open state, so it cannot approve the story.
The user requested this reference back on 2026-08-05; it was restored
byte-for-byte from commit `46ac16a` after cleanup had removed it. This is a
deliberate review reference, not a backup and not a reversal of its superseded
story status.

The accepted first Flow video remains. Future Flow generation restarts from
the second video piece with three sequential eight-second First Frame + Last
Frame clips. At the observed Lite price, the planned spend is 30 credits and
20 credits remain as retry reserve from a 50-credit refresh. Prompts stay on
hold until a v02 preview implements `../01_story/FALSIFICATION_VISUAL.md` and the complete
visual route is approved.

Current primary review:
`../renders/edit/navier-stokes-four-flow-clips-v01-32s.mp4` — 720 × 1280, 24
fps, 768 frames, exactly 32.000 seconds. It is built only from the four original
eight-second Downloads sources in story order: S01 → S02 → S05 → S06. No clip
is trimmed, retimed, slowed, frozen or given an overlay/placeholder. Source
audio is preserved. The only successful cleanup is removal of the generated
`4%/6%` label from S06; its motion, marker behavior, duration and audio remain
intact.

S05 still contains the generated amber beam/ring. Local spatial and
colour-mask removal tests visibly damaged the coral fluid, so none was used.
S05 must be regenerated without any scan, line, beam, cone or ring when credits
return. The earlier 29.70-second and 44.2-second review files are superseded;
do not use them for the current decision.

Current review file:
`../renders/edit/navier-stokes-animatic-v02-raw-plates-720x1280.mp4` — 720 ×
1280, 30 fps, 891 frames, exactly 29.700 seconds. It uses the complete raw
eight-second S02 generation across S02 and the first two seconds of S03, so its
useful ending remains visible. S05 uses a normal-speed cut from its raw source;
no synthetic slowdown or freeze is used. S01/S03-tail/S04/S06 remain explicit
still placeholders and S07 is the end card. All planned English overlays are
burned in for timing review. The audio track is intentionally silent.

The earlier `navier-stokes-animatic-v01-720x1280.mp4` was superseded because it
used locally repaired S02/S05 derivatives; its binary has been removed.

## Current shot state

| Shot | Source | Decision | Edit-safe file | Further generation |
| --- | --- | --- | --- | --- |
| S02 | `../renders/flow/drafts/S02_flow_lite_720p_take02.mp4` | Usable provisional plate with deterministic trim, baked gradient and final-frame hold; not final-locked | `../renders/flow/edits/S02_flow_lite_take02_6s_gradient_freeze.mp4` | Reassess only in the complete animatic |
| S05 | `../renders/flow/drafts/S05_flow_lite_720p_take01.mp4` | Fluid/camera are usable, but the generated amber beam and ring must not appear in the final film | Raw source is review-only | Regenerate the same macro world with no scan, line, beam, cone or ring |
| S06 | `../renders/flow/drafts/S06_flow_lite_720p_take01_rejected_text_and_markers.mp4` | Flow take 01 rejected because prompt instructions became visible text and markers morphed; deterministic still animation also rejected because the fluid froze and the synthetic zoom read as unstable | Not accepted | Keep take 02 prompt ready but do not generate before animatic review |
| S01 | `../renders/flow/drafts/S01_flow_lite_720p_take01.mp4` | Raw eight-second opening downloaded and preserved without local motion edits | Raw source itself | Inspect the entire take in the raw Flow sequence before deciding |

The rejected S02 take 01 binary has been removed. Exact acceptance fields and
defect notes remain in `REVIEW_SHEET.csv`.

## Reproducible local edits

- `../scripts/render_s02_accept.sh` reproduces the six-second S02 plate.
- `../scripts/render_s05_salvage.sh` reproduces the 3.83-second clean silent S05
  plate from the accepted clean tail of the generated source.
- `../scripts/render_s06_frontier.sh` reproduces the rejected 4.167-second silent
  S06 experiment. Keep it for audit only; do not use it in the final edit.
- S05 has no baked scan line. Its exact evidence label and pointer belong to
  the final typography pass.
- S05 source audio is omitted because slowing it exposed audible artifacts;
  final narration, ambience and effects are mixed separately.

No Kdenlive installation is required for these deterministic repairs. CapCut
Online remains optional for visual inspection, but it flattened tested
transparent PNG overlays to opaque black. Therefore any required alpha overlay
must be baked by the rendering scripts or recreated as a final vector/text
pass.

## Flow Agent behavior observed

- Keep `Confirm before generating` on `Always`.
- A prompt should initiate generation directly. Do not add `Do not generate
  yet`; that creates an unnecessary conversational confirmation.
- After the main prompt, type `Approve` as a separate Agent message to trigger
  Flow's built-in cost card. The typed word does not spend credits.
- Credits are charged only after clicking `Approve` on that built-in card.
- Never select `Approve, do not ask again`.
- The verified Lite offer is one 8-second 720p video for 10 credits.
- Four 10-credit generations have been used in this pass: S02 take 01, S02
  take 02, S05 take 01 and the rejected S06 take 01. From the observed
  50-credit daily balance, the expected remaining balance is 10 credits;
  S01 used the remaining 10 credits. The expected balance is now zero until
  the next refresh; verify the live UI before any later generation.

## Next exact action

1. Read `../01_story/FALSIFICATION_VISUAL.md` and review the v01 morph only for camera
   continuity.
2. Rebuild K06–K09 as v02 planning states: expected decay, measured persistent
   far-shell response, saved PNT-13 rejection and PNT-12 open frontier.
3. Review v02 at phone size and require every acceptance item in that file.
4. Only after approval, rewrite the Flow, OpenArt and Higgsfield prompts and
   regenerate from the second video piece. S05's amber beam/ring remains a
   rejected reference, not a local-retouch task.

Do not repair a failed Flow take with synthetic zoom, a freeze frame, geometry warp
or another local motion effect. A normal timeline cut from the eight-second
source is allowed only after the raw sequence is approved. Any camera
instability, frozen fluid or cube deformation requires a recorded rejection
and one targeted prompt revision before regeneration.

## S06 rejection lesson

Flow rendered the motion instructions `4%` and `6%` as visible labels. It also
turned several fixed checkpoint cubes into bright spike-like glyphs. The
source is therefore evidence of a failed generation, not an edit source.

Future visual prompts must avoid numeric timing, percentages and marker counts
inside the scene description. Put model, duration, resolution and output count
only in the configuration/confirmation paragraph. Describe camera movement as
`barely perceptible` and preserve `every existing marker` without asking for a
travelling marker pulse. The exact audited marker count belongs to the
canonical still and editor typography, not generative motion.

The deterministic still-animation fallback was also rejected by visual review:
it froze the fluid and the tiny synthetic zoom looked like unstable camera
motion. A valid S06 must contain clearly visible continuous fluid animation
while the camera remains level on a single centered optical axis.
