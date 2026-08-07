# Production log

This file persists the decisions made during the supervised production pass.
It is the durable record; chat instructions are not the source of truth.

Last updated: 2026-08-05  
Current next action: watch the completed caption-led post master at phone size;
do not spend the held 10-credit retry unless that review finds a defect that
cannot be repaired deterministically

## V04 review — remove the actual Flow overlap loops — 2026-08-05

V03 still exposed the persistent geometry change because a dissolve cannot
make two different states continuous. Exact frame matching found the real
structure of the downloaded Extend chain: frame 359 at 14.958 seconds recurs
almost identically at frame 383 at 15.958 seconds, and frame 551 at 22.958
seconds recurs almost identically at frame 575 at 23.958 seconds. These are
approximately one-second Flow overlap loops, not isolated bad frames.

The current review is
`../renders/edit/navier-stokes-flow-post-v04-review-720p.mp4`. It removes each
overlap from the repeated frame to its matching recurrence, then retimes the
remaining 169-frame continuation back to its eight-second timeline position.
No visual dissolve remains at 15 or 23 seconds. Contact-strip review shows the
matched states continuing without a geometry jump. The file is 720 × 1280,
30 fps, 891 frames and 29.700 seconds; full decode passed and no Flow credits
were used.

## V03 review after residual 15-second jump — 2026-08-05

User review found that V02 still appeared to jump around 14–15 seconds. The
two removed source frames were not the complete issue: the geometry on the new
side remains persistently different, so a two-output-frame blend still visibly
flipped between states. The current review is
`../renders/edit/navier-stokes-flow-post-v03-review-720p.mp4`. It keeps the same
two-frame source removal but widens only the 15-second bridge to six output
frames (14.80–15.00). Frame-strip review now shows a monotonic dissolve without
an abrupt geometry change. The accepted two-frame repair at 23 seconds remains
unchanged. V03 is 720 × 1280, 30 fps, 891 frames and 29.700 seconds; full decode
passed. No Flow credits were used.

## Native-resolution v02 review after user QC — 2026-08-05

The user rejected v01 after observing another picture jump around 13–14
seconds and audible distortion at several points. Frame-delta analysis located
the persistent source geometry change at exactly 15.000 seconds; playback and
the simultaneous event-card change can make it feel slightly earlier. V02 now
repairs both source jumps, at 15.000 and 23.000 seconds, by dropping only the
first two frames of each new state and applying a two-frame dissolve.

Current review file:
`../renders/edit/navier-stokes-flow-post-v02-review-720p.mp4` — native 720 ×
1280 picture, 30 fps, 891 frames and 29.700 seconds. It fully decodes without
error. It intentionally avoids the unnecessary 1080p upscale until picture
and sound are approved.

V01's broadband de-click processing was removed completely because it could
mistake the score's dry musical transients for defects. V02 preserves the
original score: the 7-second join receives only a 12 ms fade-out/fade-in, and
the 15- and 23-second joins use the same 83 ms equal-power crossfade as the
picture edit. No other audio repair filter is active. Measured output is about
−16.28 LUFS with −1.33 dBTP peak. If distortion remains away from those exact
joins, it is present in the Flow-generated source rather than introduced by
post.

## Completed continuous Flow chain and first post master — 2026-08-05

The supervised paid route is complete: one Veo 3.1 Lite image-to-video V01 and
three sequential Extend actions were generated in the new continuous-morph
project for 40 credits total. The fifth 10-credit slot remains reserved and was
not spent. The downloaded chain is
`/home/raxim/Downloads/Initial_Scene_-_2026-08-05_202608051411.mp4`: H.264/AAC,
720 × 1280, 24 fps and 31.000 seconds. Flow reused approximately one second of
context at each Extend, so the delivered chain is 31 rather than 32 seconds.

Raw review found continuous camera/fluid motion at the 7- and 15-second
Extend boundaries. The user identified a small geometry jump at 23 seconds.
The new post script drops only source frames 552 and 553, bridges the adjacent
states with a two-frame dissolve, and uses the generated stable tail to retain
the locked 29.700-second edit. Audio cleanup is restricted to 240 ms around
the 7- and 15-second joins; the 23-second audio follows the same short
crossfade as the picture. The rest of the generated score is unchanged apart
from a safe global gain adjustment.

Superseded caption-led post master:
`../renders/edit/navier-stokes-flow-post-v01-29.70s.mp4`. FFprobe verifies H.264
High Profile/AAC, 1080 × 1920, 30 fps, 891 frames and exactly 29.700 seconds.
The file fully decodes without errors; measured audio is about −16.35 LUFS with
−1.11 dBTP peak. It includes the immediate topic/question, sequential research
story, a separate frequency-shell diagram with no pointer into the liquid, the
expanding glass cards, one cyan-to-coral perimeter sweep and four brief final
application glints. The original Downloads file remains untouched.

The generated score is the only audio in this caption-led master. A spoken
narration layer is not baked because no recorded or selected English voice
asset exists yet. That optional voice pass must lower the score under speech
and follow `../01_story/narration/en.md`; it does not require another Flow
generation.

Preflight review, 2026-08-05: passed. The V01 anchor exists at 720 × 1280; the
planning preview is valid H.264 at 720 × 1280, 24 fps and 30.000 seconds; the
final pullback strip preserves the same coral structure as it shrinks into the
central cube. V01–V04 each contain the no-cut, no-generated-text, no-voice and
continuous-score contracts. Current-route documentation contains no active
silent-source/music-in-post contradiction, and `git diff --check` passes. The
preview remains a deterministic storyboard/timing reference, so actual camera
and music continuity must still be accepted after each paid Flow action.

Final-motion clarification, 2026-08-05: the static pullback strip demonstrates
only framing, scale and preservation of the coral structure. It is not a target
for frozen V04 scenery. The paid prompt now explicitly requires the weather
vortex to rotate, wing streamlines to travel, vessel flow to pulse, the turbine
rotor to spin and connected cyan wavefronts to propagate simultaneously while
the camera pulls back. Any static application field is rejected.

Tail-pointer correction, 2026-08-05: the preview pointer had landed on the cyan
vortex core even though no spatial point in the liquid is the mathematical
tail or a raw C++ field location. Moving it to another decorative feature would
still imply false precision, so the pointer, circle and fluid target were
removed completely. The tail now appears only as the far end of the separate
frequency-shell diagram, explicitly labelled as not a physical wave tail.

Far-shell transition correction, 2026-08-05: the earlier no-credit preview
crossfaded from the legacy K05 still to a differently framed K06 base, making
the cube appear to jump backward when the diagram arrived. K04–K09 now use
progressively tighter crops of one source image. The diagram changes as an
editor overlay while the planned camera continues forward; the paid Flow route
uses Extend on the same moving video and must not reproduce a still-image
crossfade.

Overlay-motion decision, 2026-08-05: major caption cards now use an expanding
dark-glass box, one restrained border light sweep and a rapid top-to-bottom
full-line reveal. The final pullback adds brief tracked glints to the moving
weather vortex, wing airflow, vessel flow and turbine. This is deterministic
editor animation in `../04_editorial/OVERLAY_MOTION_PROMPT.md`, not part of the
Flow prompt, because Veo must not render the exact text or UI.

Opening-title decision, 2026-08-05: the finished edit must identify the topic
immediately over the calm V01 opening with `NAVIER–STOKES REGULARITY` and the
question `CAN A SMOOTH FLUID BECOME INFINITELY ROUGH?`. This is deterministic
post typography and does not justify regenerating an otherwise accepted V01.

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

The earlier S01 remains as an audit/reference take, but the new paid route
starts fresh so all four slots belong to one extension chain. Generate one
eight-second Veo 3.1 Lite V01, then use Extend on that accepted video three
times for V02–V04. The four actions are fixed slots, not scenes. At the current
non-Ultra price the planned spend is 40 credits; the remaining 10 credits from
a 50-credit refresh are reserved for one targeted retry.

The user reported a refreshed credit balance on 2026-08-05. Use a new project
named `Navier–Stokes — Continuous Morph V01–V04 — 9x16`; do not approve the
stale 10-credit card still visible in the legacy `Lemma Stress Test` project.
The live confirmation card must still show one Lite generation at no more than
10 credits before each paid action.

Audio decision, 2026-08-05: the earlier independent-shot prompts intentionally
requested very soft glassy-water ambience or ambient electronic pulses, but
they did not define a single full-film musical score. For the current Extend
chain, V01 now generates one restrained premium science-documentary electronic
cue and V02–V04 must continue that exact cue without a restart, silence gap,
tempo/key shift or instrument change. No separate music track is planned. The
generated footage must contain no voice, speech, singing or lyrics. A broken
audio seam is a rejection condition and can justify the held 10-credit retry.
If Flow exposes a generated-audio toggle, it must remain enabled throughout
the chain; reject a silent V01. An explicit `Audio Generation Failed` error is
expected by Google's current help to refund credits, but verify the live
balance before retrying.

Current no-credit story preview:
`../renders/edit/navier-stokes-10-frame-story-previs-v02-30s.mp4` — 720 ×
1280, 24 fps, 720 frames and exactly 30.000 seconds. K06–K09 move progressively
deeper into the same fluid while editor-owned overlays distinguish expected
decay, the persistent far-shell response, the isolated obstruction, the saved
PNT-13 rejection and the PNT-12 open state. The final three seconds use a true
scale pullback on the same planning image: the same coral structure remains
visible in the central cube while weather, flight, blood-flow and energy forms
are revealed around it. This is still planning footage, not a generated Flow
take. User approval is required before paid generation.

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

1. Review
   `../renders/edit/navier-stokes-10-frame-story-previs-v02-30s.mp4` at phone
   size against every item in `../01_story/FALSIFICATION_VISUAL.md`.
2. Only after user approval, follow
   `../03_platforms/flow/CHAIN_32S_RUNBOOK.md`: V01 plus three Extend actions,
   40 credits total and a 10-credit retry reserve. S05's amber beam/ring remains
   a rejected reference, not a local-retouch task.

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
