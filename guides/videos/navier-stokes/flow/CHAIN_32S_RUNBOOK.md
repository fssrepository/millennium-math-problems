# Flow — continuous 32-second chain from three replacement generations

> **HOLD — DO NOT PASTE THESE PROMPTS YET.** Review the ten-frame preview in
> `../renders/edit/navier-stokes-10-frame-morph-previs-v01-30s.mp4` only for
> camera continuity, then read `../FALSIFICATION_VISUAL.md`. Rewrite the
> prompts and endpoint choices only after a v02 preview visibly includes the
> PNT-13 obstruction and PNT-12 state handoff.

Keep the accepted eight-second first Flow video. Restart generation from the
second video piece with three eight-second **First Frame + Last Frame** clips.
The final chain remains exactly 32 seconds.

| Timeline | Source | First frame | Last frame |
| --- | --- | --- | --- |
| 00:00–00:08 | accepted existing S01 | existing generated motion | actual S01 final frame |
| 00:08–00:16 | C01 replacement | actual S01 final frame | endpoint to lock after previs review |
| 00:16–00:24 | C02 new transition | C01 actual saved final frame | `chain-anchor-02-s05-last-clean.png` target |
| 00:24–00:32 | C03 new transition | C02 actual saved final frame | impact endpoint to lock after previs review |

The accepted original S01 remains. The defective S05 scan/beam and original
S06 percentage counter do not enter the replacement chain. The visual world
stays continuous. Only the
liquid and camera may move; no new diagram, interface or symbolic event may
appear.

The global camera-energy curve is mandatory: existing S01 is fastest, C01 is
slower, C02 is slower again and C03 has the lowest camera energy. No later clip may
accelerate above the previous clip. The final height-field reveal gains visual
impact from parallax and depth, not from late speed.

## Boundary-continuity rule

The three replacement generations are sequential, not independent. Save the
accepted S01's exact final frame and use it as C01's First Frame. After C01 and
C02, always use the previous accepted clip's exact saved final frame as the
next First Frame.

Do not start a later clip from a pre-extracted target merely because it looks
similar. The actual output frame is the continuity handoff. This guarantees
matching pixels at the clip
boundary, although a generative model cannot guarantee correct motion inside a
clip. Inspect every result for frozen liquid, a speed discontinuity or a camera
direction change before spending on the next clip.

## Credit and configuration gate

Run exactly three generations after the previs and prompts are approved. For
each one require:

- media type: exactly one video;
- input mode: Frames to Video, with both **First Frame** and **Last Frame**;
- model: Veo 3.1 Lite;
- portrait 9:16;
- 720p;
- duration: 8 seconds;
- outputs: one;
- displayed cost: no more than 10 credits.

Reject the proposal if either image is described only as an ingredient or
reference. Reject a model switch, missing end frame, multiple outputs or a
cost above 10 credits.

The expanded Agent uses two approvals:

1. attach the two images and paste the complete prompt;
2. after the Agent describes the job, type `Approve` as a separate message;
3. inspect the built-in cost card;
4. click the first built-in `Approve` only when every field above matches.

Never click `Approve, do not ask again`. At the observed price, three accepted
runs cost 30 credits and leave 20 credits from a 50-credit refresh as reserve.

## C01 — cube to living-fluid macro

Attach these exact files in these roles:

- First Frame: `../assets/flow-chain/chain-anchor-00-s01-last.png`
- Last Frame: `../assets/flow-chain/chain-anchor-01-s02-last.png`

Paste this complete prompt:

```text
Create exactly one continuous video. Use the first attached image as the exact
first frame and the second attached image as the exact last frame. They are
timeline endpoints, not ingredients or loose visual references. Do not create
or edit a still image.

Use Veo 3.1 Lite, portrait 9:16, 720p, 8 seconds and one output. Prepare this
exact generation, then wait for the user's separate single-word message
`Approve`. When it arrives, invoke the video-generation action so Flow shows
its mandatory built-in credit confirmation.

Create one uninterrupted optical push-in from the full glass cube in the first
frame to the interior fluid macro in the last frame. The camera travels
smoothly through the transparent cube boundary without impact, splash or
surface disturbance, then follows the same coral-red fluorescent liquid dye
inside the same coherent teal liquid. Continue the apparent incoming camera
movement immediately with no opening pause. Acceleration is confident at first
and eases into a steady, slow forward camera velocity through the exact final
macro composition; do not stop or hold on the last frame. Preserve physically
plausible continuous liquid motion, the identity and color of the coral dye,
the teal flow, the rigid sealed cube and the dark premium science-documentary
world. Only the camera and existing liquid move. The first and last
compositions must match the attached frames exactly. No cut and no voice.

Do not create a scan, line, vertical line, horizontal line, beam, plane, cone,
wedge, arrow, flashlight, ring, circle, aperture, pulse, progress indicator or
moving overlay. No visible text or symbols. No letters, numbers, percentage
signs, percentage counters, equations, subtitles, logos, watermark, HUD or UI.
No people, rope, braid, solid coral object, extra coral filaments, duplicated
or deforming cube, broken glass, splash, explosion, outward particles, camera
shake, micro-wobble, flicker, random cut, geometry morphing, color drift or
frozen liquid.
```

Save as `C01_flow_lite_first-last_8s_take01.mp4`.

After C01 passes review, open C01, pause on its actual final frame, choose
**Save frame**, and use that newly saved frame for C02. This is mandatory.

## C02 — continuous macro passage, no scan effect

Attach these exact files in these roles:

- First Frame: C01's actual final frame saved from the generated video
- Last Frame: `../assets/flow-chain/chain-anchor-02-s05-last-clean.png`

Paste this complete prompt:

```text
Create exactly one continuous video. Use the first attached image as the exact
first frame and the second attached image as the exact last frame. They are
timeline endpoints, not ingredients or loose visual references. Do not create
or edit a still image.

Use Veo 3.1 Lite, portrait 9:16, 720p, 8 seconds and one output. Prepare this
exact generation, then wait for the user's separate single-word message
`Approve`. When it arrives, invoke the video-generation action so Flow shows
its mandatory built-in credit confirmation.

Remain inside one continuous macro fluid world. Track smoothly along the same
single coral-red fluorescent liquid-dye structure as it stretches, folds and
relaxes inside the same coherent teal liquid. For the opening second, continue
the incoming forward camera direction and speed without a pause. Then use a
controlled macro arc with natural parallax. During the final two seconds,
transition into a gentle backward camera drift and carry that same backward
velocity through the exact last frame; do not stop or hold. Preserve fluid
continuity, soft diffusing dye
edges, physically plausible liquid motion, the dark cyan-and-coral palette and
premium science-documentary realism. Only the camera and existing liquid move.
The first and last compositions must match the attached frames exactly. No cut
and no voice.

Do not create a scan, line, vertical line, horizontal line, beam, plane, cone,
wedge, arrow, flashlight, ring, circle, aperture, pulse, progress indicator or
moving overlay. No visible text or symbols. No letters, numbers, percentage
signs, percentage counters, equations, subtitles, logos, watermark, HUD or UI.
No people, rope, braid, solid coral object, extra coral filaments, glass break,
splash, explosion, checkmark, victory image, camera shake, micro-wobble,
flicker, random cut, geometry morphing, color drift or frozen liquid.
```

Save as `C02_flow_lite_first-last_8s_take01.mp4`.

After C02 passes review, open C02, pause on its actual final frame, choose
**Save frame**, and use that newly saved frame for C03. This is mandatory.

## C03 — macro to frontier reveal

Attach these exact files in these roles:

- First Frame: C02's actual final frame saved from the generated video
- Last Frame: `../assets/keyframe-06-frontier.png`

Paste this complete prompt:

```text
Create exactly one continuous video. Use the first attached image as the exact
first frame and the second attached image as the exact last frame. They are
timeline endpoints, not ingredients or loose visual references. Do not create
or edit a still image.

Use Veo 3.1 Lite, portrait 9:16, 720p, 8 seconds and one output. Prepare this
exact generation, then wait for the user's separate single-word message
`Approve`. When it arrives, invoke the video-generation action so Flow shows
its mandatory built-in credit confirmation.

Continue the exact incoming backward camera direction and speed from the first
frame with no pause. Perform one visually dramatic but smooth macro-to-wide
optical pullback. Begin inside the
same living coral-red fluorescent liquid dye and coherent teal liquid from the
exact first frame. Pull backward through the rigid transparent cube boundary
without impact or disturbance, then continue outward until the central cube,
nested transparent cyan shells and restrained curved amber evidence path are
revealed in the exact final composition. The same liquid remains visibly alive
and identifiable throughout the move. Reveal existing endpoint structures by
camera motion and occlusion only; do not pop, spawn, count or transform them.
Only the camera and existing liquid move. Use strong readable parallax through
the shell layers, then decelerate only during the final second into a perfectly
level, centered, stable final frame. Preserve the dark premium science-documentary
world. No cut and no voice.

Do not create a scan, line, vertical line, horizontal line, beam, plane, cone,
wedge, arrow, flashlight, ring, circle, aperture, pulse, progress indicator or
moving overlay. No visible text or symbols. No letters, numbers, percentage
signs, percentage counters, equations, subtitles, logos, watermark, HUD or UI.
No people, extra markers, travelling marker pulse, rope, braid, solid coral
object, extra coral filaments, duplicated or deforming cube, broken glass,
splash, explosion, checkmark, victory flare, free-floating sparks, camera
shake, micro-wobble, flicker, random cut, geometry morphing, color drift or
frozen liquid.
```

Save as `C03_flow_lite_first-last_8s_take01.mp4`.

## After each download

Watch the result once at normal speed and once at 0.5×. Reject it if either
endpoint is wrong, a cut appears, fluid freezes, the cube deforms, a scan line
appears, or text/numbers/percentages are generated. Do not repair a failed
transition with synthetic zoom, slowdown, freeze or geometry warp.

Place accepted downloads in `../renders/flow/chain/`, then run:

```bash
./guides/videos/navier-stokes/scripts/build_flow_chain_32s.sh
```

The review output will be
`../renders/edit/navier-stokes-flow-chain-v01-32s.mp4`.
