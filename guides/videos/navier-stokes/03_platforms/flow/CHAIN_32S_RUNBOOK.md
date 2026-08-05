# Flow — one continuous 32-second camera morph in four Lite slots

> **READY FOR FINAL PREVIEW REVIEW — DO NOT SPEND YET.** Watch
> `../../renders/edit/navier-stokes-10-frame-story-previs-v02-30s.mp4` at phone
> size. Generate only after the user confirms that the camera path and final
> application morph read as one uninterrupted video.

Create a clean Flow project named
`Navier–Stokes — Continuous Morph V01–V04 — 9x16`. Do not reuse the older
`Lemma Stress Test` project: it contains legacy independent takes and stale
approval cards. Flow credits are account-wide, so the new project does not
create or consume a separate balance.

## What “continuous” means

This is not a four-scene storyboard. It is one Veo-generated clip extended
three times:

| Slot | Flow action | Running duration | Camera beat | Cost ceiling |
| --- | --- | ---: | --- | ---: |
| V01 | Generate one 8 s clip | 8 s | outside cube → push toward fluid | 10 |
| V02 | Extend V01 by 8 s | 16 s | accelerate inside the same liquid | 10 |
| V03 | Extend V02 by 8 s | 24 s | decelerate onto the observed tail proxy | 10 |
| V04 | Extend V03 by 8 s | 32 s | near-arrest → smooth scale pullback and application reveal | 10 |

The planned spend is 40 credits. From a 50-credit daily balance, keep the last
10 credits untouched as one targeted retry. Do not request variants or multiple
outputs: Flow charges per generation, not per written request.

Google's current Flow documentation says Veo 3.1 Lite supports 8-second videos
and Extend, and lists a 10-credit cost per generation for non-Ultra users.
Costs can change; the live confirmation card is authoritative:

- https://support.google.com/flow/answer/16352836
- https://support.google.com/flow/answer/16526234
- https://support.google.com/flow/answer/16935718

## Non-negotiable continuity lock

- Use **Veo 3.1 Lite**, portrait **9:16**, **720p**, **8 seconds**, **one
  output** for all four paid actions.
- V02, V03 and V04 must use **Extend** on the current accepted video. Never
  create a fresh clip from a similar still, never use Scene Builder and never
  ask for four separate scenes.
- One camera has one uninterrupted trajectory. It may accelerate, decelerate
  and zoom. After the tail inspection it may pass smoothly through zero
  velocity into one deliberate pullback, but it never cuts, teleports or resets
  its optical axis at a slot boundary.
- The exact same coral-red liquid dye and teal fluid persist through every
  extension. Their motion never freezes at a slot boundary.
- During the late pullback, the same central coral response remains visible and
  shrinks into the central cube while weather, wing, vessel and turbine flow
  forms are revealed in the widening field. They do not replace it or pop in
  as montage panels.
- Generate no text, numbers, percentages, equations, labels, pointer, ledger,
  UI, subtitles or narration. Those are deterministic editor overlays.
- Generate one continuous original instrumental score together with V01 and
  preserve that exact score through V02–V04. Its identity is restrained
  premium science-documentary electronica: a steady mid-tempo pulse, warm low
  synth, delicate glassy-water texture, sparse dry ticks and a soft airy
  whoosh. It is intelligent and curious, tense but not ominous.
- At every Extend boundary, keep the same tempo, tonal center, instrument
  palette, ambience and musical phrase. Never restart the cue, introduce a new
  downbeat or fade to silence between slots. No dialogue, narration, spoken
  words, singing, lyrics or choir.
- No heroic trailer rise, orchestral boom, EDM drop, pop hook or triumphant
  cadence. The score follows the camera's acceleration, near-arrest and final
  pullback without becoming a separate musical sequence.

## Credit gate before every action

Proceed only when Flow's built-in card shows all of the following:

- Veo 3.1 Lite;
- portrait 9:16;
- 720p;
- 8 seconds;
- one output;
- 10 credits or less.

Reject the card if it proposes multiple outputs, Fast, Quality, Omni Flash, a
fresh scene, a still image, an ingredient-based composition or more than 10
credits. Never choose `Approve, do not ask again`.

After each slot, watch the boundary at normal speed and 0.5× before spending
the next 10 credits. Reject on a cut, frozen liquid, camera reset, direction
change, duplicated coral dye, geometry morph before the final application beat,
visible text or any new interface graphic. Listen on headphones as well: reject
an audio click, silence gap, cue restart, tempo or key jump, changed instrument
palette, sudden loudness step, voice, singing or lyrics. An otherwise good
visual with a broken musical seam is the one justified use of the held retry.

## V01 — initial generation

Configuration:

- mode: Frames to Video → First Frame only;
- First Frame: `../../assets/flow-chain/chain-anchor-start-s01-first.png`;
- no Last Frame and no Ingredients.

Paste this prompt:

```text
Create the opening eight seconds of one uninterrupted portrait science-
documentary camera move. This is the beginning of a single video that will be
extended three times, not a standalone scene.

Begin exactly from the uploaded first frame. Preserve the rigid transparent
glass cube, the same coherent luminous teal incompressible liquid, and the one
recognizable coral-red fluorescent liquid-dye vortex with soft diffusing edges
and living internal wisps. The near-black research space, cyan rim light and
subtle floor reflections remain physically stable.

The camera begins calm and level, then pushes forward toward the cube with
steadily increasing confidence. Approach the glass on one centered optical
axis and begin entering the liquid without impact, splash, broken glass or a
pause. Carry clear forward velocity through the final frame so an extension
can continue the same move immediately. The teal and coral liquid remain alive
throughout. Only the camera and the existing liquid move. No cut.

Create one original continuous instrumental score as part of the generated
audio. Use restrained premium science-documentary electronica: a steady
mid-tempo pulse, warm low synth, delicate glassy-water textures, sparse dry
percussive ticks and a soft airy whoosh. Keep it intelligent and curious,
slightly tense but not ominous. Begin sparsely and let its energy rise gently
with the forward camera motion. This is the opening of one score that the next
three Extends must continue, so carry the musical phrase and ambience through
the final audio sample without a cadence, fade or hard stop.

No visible text or symbols. No letters, numbers, percentages, equations,
subtitles, logos, watermark, HUD or UI. No people, extra objects, scan line,
beam, ring, radiating visual pulse, pointer, ledger, markers, scene transition,
portal, camera shake, micro-wobble, flicker, frozen liquid, duplicated dye,
deforming cube, geometry jump or color drift. No dialogue, narration, spoken
words, singing, lyrics, choir, heroic trailer rise, orchestral boom, EDM drop
or pop hook.
```

Save the accepted result as `V01_flow_lite_continuous_8s_take01.mp4`.

## V02 — first Extend

Select the accepted V01 video and click **Extend**. Do not upload a new start
frame. Paste this continuation prompt:

```text
Continue the existing video with no cut, pause or camera reset. Inherit the
exact final-frame pixels, forward camera direction, velocity, optical axis,
lighting, coral-red dye identity and teal liquid motion from the current clip.

The same camera passes fully inside the same liquid and accelerates along the
same coral dye structure. Use a controlled macro flight with strong natural
parallax: follow its folds and living wisps while coherent teal streams move
around it. The move becomes the most energetic part of the film, then eases
slightly without stopping. Preserve one physical fluid simulation and carry
forward velocity through the final frame so the next Extend continues it.
Only the camera and existing liquid move. No cut and no new location.

Continue the exact existing instrumental score seamlessly from the current
last audio sample. Preserve its tempo, tonal center, pulse, warm low synth,
glassy-water texture, dry ticks, airy whoosh and acoustic space. Increase its
rhythmic density with the camera acceleration, then ease slightly without a
new intro, downbeat, cue change, cadence, fade, silence gap or loudness jump.

No visible text or symbols. No letters, numbers, percentages, equations,
subtitles, logos, watermark, HUD or UI. No people, extra objects, scan line,
beam, ring, radiating visual pulse, pointer, ledger, markers, scene transition,
portal, camera shake, micro-wobble, flicker, frozen liquid, duplicated dye,
geometry jump or color drift. No dialogue, narration, spoken words, singing,
lyrics, choir, heroic trailer rise, orchestral boom, EDM drop or pop hook.
```

Save the accepted extended result as
`V02_flow_lite_continuous_16s_take01.mp4`.

## V03 — second Extend

Select the accepted V02 extended video and click **Extend**. Paste:

```text
Continue the same uninterrupted camera move and fluid simulation directly from
the current final frame. Do not restart the shot. Keep the exact same
coral-red liquid dye, teal flow, dark environment, camera axis and forward
direction.

Follow the same coral structure deeper into macro detail. Gradually decelerate
from the earlier fast flight into a precise close inspection of one identifiable
coral fold that will serve as the editor's visual proxy for the far-shell
response. The fluid stays visibly alive while the camera becomes slow and
deliberate. Use shallow parallax and a gentle optical push-in; never freeze,
pull back or hold. Carry a small but unmistakable forward velocity through the
final frame. Only the camera and existing liquid move. No cut.

Continue the exact same score and ambience with no musical boundary. Gradually
thin the dry ticks and pulse as the camera slows, leaving the warm low tone,
glassy liquid texture and airy motion nearly exposed at the close inspection.
Approach near-silence without ever stopping the audio bed. Keep the same tempo,
tonal center and instruments; no new cue, intro, downbeat, cadence, silence gap
or loudness jump.

Do not generate the analytical guide or its labels: no spectral bands, graph,
pointer, status marker, text, numbers, percentages, equations, subtitles,
logos, watermark, HUD or UI. No scan, beam, ring, radiating visual pulse, scene
transition, portal, new location, camera shake, micro-wobble, flicker, frozen
liquid, duplicated dye, geometry jump or color drift. No dialogue, narration,
spoken words, singing, lyrics, choir, heroic trailer rise, orchestral boom, EDM
drop or pop hook.
```

Save the accepted extended result as
`V03_flow_lite_continuous_24s_take01.mp4`.

## V04 — third Extend and continuous application morph

Select the accepted V03 extended video and click **Extend**. Paste:

```text
Continue directly from the current final frame as the same uninterrupted video.
Inherit the exact coral feature, teal flow, forward camera direction, optical
axis, lighting and living fluid velocity. No cut, pause or new scene.

Complete the close inspection while the camera still creeps forward. Slow
smoothly through zero velocity, then perform one fast but controlled optical
pullback on the same centered axis. Keep the exact coral feature visible and
identifiable throughout: it becomes smaller inside the same central
transparent cube as the field of view expands. Never replace it with a new
object.

The widening view reveals one connected flow world already growing from the
same cyan and coral streamlines: a weather-system spiral above, airflow
wrapping an aircraft wing, branching blood-flow geometry, and a turbine flow
below. The applications arrive through camera scale, parallax and continuous
fluid transformation around the preserved central cube. This is one morph
observed by one moving camera, not a montage, grid, sequence of scenes or set
of portals. Ease only near the final frame while every flow remains alive.

Continue the exact existing score from its current audio sample with no seam.
Bring back a little of the same restrained pulse and widen the same harmonic
space as the camera pulls back, while preserving the tonal center, instrument
palette, glassy-water texture and airy whoosh. During the final portion settle
into a long stable sustained tail that remains clean if the release edit trims
the master slightly. Do not start a new cue, restart the rhythm, change genre,
add a trailer rise or create a triumphant final cadence.

No visible text or symbols. No letters, numbers, percentages, equations,
subtitles, logos, watermark, HUD or UI. No people, panel borders, split screen,
hard cut, cross-cut, scene reset, portal, object pop-in, explosion, camera
shake, micro-wobble, flicker, freeze, duplicated subjects, color drift or
voice. No dialogue, narration, spoken words, singing, lyrics or choir.
```

Save the accepted final extended result as
`V04_flow_lite_continuous_32s_take01.mp4`.

## Reserve rule

Do not spend the fifth 10-credit slot automatically. Use it only when one
specific defect makes the 32-second result unusable. Retry the failed slot or
extension once with one targeted wording change; do not regenerate all four.
If the same defect remains, stop and review rather than consuming more credit.

## Final review

Review the final 32-second result three times:

1. with no overlays, to judge camera and fluid continuity;
2. on headphones, to judge the single continuous score and all three Extend
   boundaries;
3. with the planned editorial overlay timing, to judge the research story.

The footage passes only when all four slot boundaries are invisible, the tail
inspection is readable, the final applications arrive by morph rather than by
scene change, the score behaves as one uninterrupted cue, and no
generator-made text or symbolic graphics appear.

Keep the complete 32-second master. The release edit may trim only the tail to
the locked 29.70-second duration; never cut around a failed slot boundary to
hide a continuity defect.
