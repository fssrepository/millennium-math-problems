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
| V04 | Extend V03 by 8 s | 32 s | resume motion → continuous application morph | 10 |

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
- One camera has one uninterrupted forward trajectory. It may accelerate,
  decelerate and zoom, but it never cuts, teleports, resets its optical axis or
  reverses merely to start a new slot.
- The exact same coral-red liquid dye and teal fluid persist through every
  extension. Their motion never freezes at a slot boundary.
- The late weather, wing, vessel and turbine forms grow by a continuous fluid
  morph from the existing streamlines. They do not pop in as montage panels.
- Generate no text, numbers, percentages, equations, labels, pointer, ledger,
  UI, subtitles or narration. Those are deterministic editor overlays.
- Prefer a silent source. At minimum: no dialogue, no voice and no music.

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
visible text or any new interface graphic.

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

No visible text or symbols. No letters, numbers, percentages, equations,
subtitles, logos, watermark, HUD or UI. No people, extra objects, scan line,
beam, ring, pulse, pointer, ledger, markers, scene transition, portal, camera
shake, micro-wobble, flicker, frozen liquid, duplicated dye, deforming cube,
geometry jump, color drift or voice.
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

No visible text or symbols. No letters, numbers, percentages, equations,
subtitles, logos, watermark, HUD or UI. No people, extra objects, scan line,
beam, ring, pulse, pointer, ledger, markers, scene transition, portal, camera
shake, micro-wobble, flicker, frozen liquid, duplicated dye, geometry jump,
color drift or voice.
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

Do not generate the analytical guide or its labels: no spectral bands, graph,
pointer, status marker, text, numbers, percentages, equations, subtitles,
logos, watermark, HUD or UI. No scan, beam, ring, pulse, scene transition,
portal, new location, camera shake, micro-wobble, flicker, frozen liquid,
duplicated dye, geometry jump, color drift or voice.
```

Save the accepted extended result as
`V03_flow_lite_continuous_24s_take01.mp4`.

## V04 — third Extend and continuous application morph

Select the accepted V03 extended video and click **Extend**. Paste:

```text
Continue directly from the current final frame as the same uninterrupted video.
Inherit the exact coral feature, teal flow, forward camera direction, optical
axis, lighting and living fluid velocity. No cut, pause or new scene.

Complete the close inspection while the camera still creeps forward, then
accelerate smoothly through the coral feature. As the camera gains speed, let
the existing cyan and coral streamlines transform continuously and physically
into one connected wide composition: a weather-system spiral, airflow wrapping
an aircraft wing, branching blood-flow geometry and a turbine flow. Each form
must grow from the same moving streamlines and remain connected in one dark
scientific world. This is one fluid morph observed by one moving camera, not a
montage, grid, sequence of scenes or set of portals. Use camera momentum,
parallax and scale change to create the spectacle. Ease only near the final
frame while the flow remains alive.

No visible text or symbols. No letters, numbers, percentages, equations,
subtitles, logos, watermark, HUD or UI. No people, panel borders, split screen,
hard cut, cross-cut, scene reset, portal, object pop-in, explosion, camera
shake, micro-wobble, flicker, freeze, duplicated subjects, color drift or
voice.
```

Save the accepted final extended result as
`V04_flow_lite_continuous_32s_take01.mp4`.

## Reserve rule

Do not spend the fifth 10-credit slot automatically. Use it only when one
specific defect makes the 32-second result unusable. Retry the failed slot or
extension once with one targeted wording change; do not regenerate all four.
If the same defect remains, stop and review rather than consuming more credit.

## Final review

Review the final 32-second result twice:

1. with no overlays, to judge camera and fluid continuity;
2. with the planned editorial overlay timing, to judge the research story.

The footage passes only when all four slot boundaries are invisible, the tail
inspection is readable, the final applications arrive by morph rather than by
scene change, and no generator-made text or symbolic graphics appear.

Keep the complete 32-second master. The release edit may trim only the tail to
the locked 29.70-second duration; never cut around a failed slot boundary to
hide a continuity defect.
