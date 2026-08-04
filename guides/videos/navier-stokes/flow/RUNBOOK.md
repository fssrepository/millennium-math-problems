# Google Flow — cost-efficient button-by-button workflow

## Current primary route

Do not spend the next credit refresh yet. First review the ten-frame,
30-second visual plan at
`../renders/edit/navier-stokes-10-frame-morph-previs-v01-30s.mp4`. Keep the
accepted first video, then rewrite and generate only the second, third and
fourth eight-second pieces from `CHAIN_32S_RUNBOOK.md`. The current chain
prompts are explicitly on hold until the previs is approved.

The standard Flow path is `Video → Frames → First Frame + Last Frame`. If the
expanded Agent is the only available interface, attach both images and approve
only when its built-in card explicitly assigns them as the first and last
frames. Two generic references are not equivalent.

This is the primary route. It uses only Flow's standard prompt box and built-in
**Frames** image-to-video mode. Generate one 720p result, inspect it, and spend
again only when the result identifies a concrete next action.

If Flow opens an expanded Agent session and the standard prompt box is not
available, use `AGENT_READY_TO_PASTE.md` for the current supervised scout. Its
exact combined prompt explicitly requests one first-frame video and preserves
the `Always confirm` credit gate. Do not place the image in `Scenes`: that area
is for assembling generated clips. Dropping a file onto the project canvas may
only import it into `All Media`; attach it to the Agent with
`+ → select asset → Add to Prompt`.

Keep Agent settings on `Confirm before generating → Always`. After sending the
complete shot prompt, wait for the Agent to describe the prepared job, then
type exactly `Approve` as a separate prompt message. This typed word triggers
the media action but does not itself spend credits. Inspect the subsequent
built-in cost card, then click its `Approve` option only when model, format,
output count and credits match the runbook. Never choose `Approve, do not ask
again`. Do not add `Do not generate yet`; that creates an unnecessary extra
confirmation exchange. If an older conversation shows `Proceed with
generation`, select it before performing the same built-in cost-card check.

On the account verified on 2026-08-04, the expanded Agent offered the S02 scout
as either Veo 3.1 Lite at 8 seconds for 10 credits or Omni Flash at 6 seconds
for 10 credits. Choose Lite at 8 seconds, instruct it to finish the camera move
by second 6 and hold the endpoint for 2 seconds, then trim the source to the
6-second S02 edit. The standard prompt-box feature matrix may expose different
durations; always obey the displayed confirmation rather than assuming parity
between the two interfaces.

The Omni Flash option is a **model switch**, not merely a way to select six
seconds. Never change model and duration silently. If a model switch becomes
necessary, give it a separate take filename and record the old model, new
model, reason, generated duration and credit charge in the review-sheet notes.
Do not mix a switched-model result into a same-model retry comparison.

### Provisional S02 animatic source

The accepted take 02 source is
`../renders/flow/drafts/S02_flow_lite_720p_take02.mp4`. The edit-safe derivative
is
`../renders/flow/edits/S02_flow_lite_take02_6s_gradient_freeze.mp4`: it bakes
the subtle upper-screen gradient into the image, keeps source motion through
5.375 seconds, freezes the last approved frame for 0.625 seconds, fades the
source audio, and ends at exactly 6.000 seconds. Import this single MP4 into
CapCut Online; do not add the transparent gradient PNG as another track.

This baked derivative became necessary because CapCut Online displayed both
the indexed-alpha and true-RGBA PNG tests as opaque black cards on 2026-08-04.
The generated source remains preserved so a different editor can reproduce or
revise the composite without another Flow generation.

### Provisional S05 animatic salvage

S05 take 01 has a stable fluid plate, but its generated scan mutates between
approximately 1.8 and 5.8 seconds into a triangular beam with a circular
handle. Preserve the source as
`../renders/flow/drafts/S05_flow_lite_720p_take01.mp4`, but do not use that
malformed scan and do not spend another Flow generation on the same plate.

Run `../scripts/render_s05_salvage.sh`. It slows the clean 6.5--8.0 second tail
to the final 3.83-second S05 duration. The ready base plate is
`../renders/flow/edits/S05_flow_lite_take01_3p83s_clean_silent.mp4`. No
replacement scan is baked into this plate: the exact evidence pointer belongs
to the later typography pass. Generated source audio is omitted because
time-stretching exposed audible artifacts; final audio is mixed later.

### Rejected S06 deterministic experiment

The Flow take is preserved as
`../renders/flow/drafts/S06_flow_lite_720p_take01_rejected_text_and_markers.mp4`
for audit only. Flow rendered numeric prompt instructions as visible labels and
morphed fixed checkpoint markers, so it is not an edit source.

Run `../scripts/render_s06_frontier.sh` only to reproduce the rejected local
experiment at
`../renders/flow/edits/S06_keyframe_4p16s_dolly_silent.mp4`. It animates the
canonical frontier keyframe with a synthetic centered dolly-out, but it freezes
the fluid and the small zoom reads as unstable camera movement. Do not use it
in the final edit.

Keep the simplified S06 take 02 prompt in `AGENT_READY_TO_PASTE.md` ready, but
do not spend on it until the complete animatic identifies S06 as the
highest-impact retry.

## Do not open

- Agent when the standard prompt box is available; the documented expanded
  Agent fallback is allowed when it is the only route shown;
- Storyboard Studio or a scene builder;
- Ingredients;
- character, avatar, or voice references;
- Extend or automatic story generation.

## Initial setup

1. Open https://flow.google and sign in.
2. Create a new project.
3. Name it `Navier Stokes — Lemma Stress Test — 9x16`.
4. Locate the large prompt box near the center or bottom of the page.
5. If **Agent** is enabled, turn it off.
6. Click the model name beside the prompt box.
7. Select **Video**.
8. Select **Frames** as the input type.
9. Drag the shot's PNG onto `+ Add start frame`.
10. Leave `+ Add end frame` empty.
11. Select:
    - aspect ratio: `9:16`;
    - model: `Veo 3.1 Lite`;
    - resolution: `720p` when a resolution control is shown;
    - outputs: `1`;
    - generation length: the nearest available value listed below.
12. Choose ambient/no speech for audio, or turn sound off.

If a label changes slightly, the path is still: **standard prompt → Video →
Frames → Add start frame**. If Lite does not expose start-frame input in your
account, use Fast with the same one-output rule.

## Scout in this order

The order minimizes wasted generations by testing the hardest camera behavior
first. It is not the final edit order. The current pass is paused for a full
animatic review: S02 and S05 are provisional plates, while S06 has no accepted
motion plate.

| Scout order | Shot | Start frame | Length | Prompt to paste |
| ---: | --- | --- | ---: | --- |
| 1 | S02 | current `../assets/keyframe-01-world.png` | 6 s used from an 8 s Agent source | `AGENT_READY_TO_PASTE.md` → S02 take 02 block |
| 2 | S05 | `../assets/keyframe-03-checkpoint.png` | 8 s | `READY_TO_PASTE.md` → complete S05 block |
| 3 | S06 | `../assets/keyframe-06-frontier.png` | 8 s | `READY_TO_PASTE.md` → complete S06 block |
| 4 | S01 | `../assets/keyframe-01-world.png` | 5 s | `READY_TO_PASTE.md` → complete S01 block |
| 5 | S03 | `../assets/keyframe-03-checkpoint.png` | 5 s | `READY_TO_PASTE.md` → complete S03 block |
| 6 | S04 | `../assets/keyframe-01-world.png` | 8 s | `READY_TO_PASTE.md` → complete S04 block |

## Exact process for one shot

Example: S02.

The original take 01 is preserved as
`../renders/flow/drafts/S02_flow_lite_720p_take01_rejected.mp4`. It must not be
used in the edit. For the targeted retry:

1. Upload the current repository copy of `keyframe-01-world.png` as a new Flow
   asset, then attach that new upload with
   `+ → select asset → Add to Prompt` in the expanded Agent route. Do not use
   Flow's cached earlier upload with the same filename.
2. Open `AGENT_READY_TO_PASTE.md` and find `S02 take 02`.
3. Copy the complete block between the three backticks.
4. Paste it into Flow's main prompt box.
5. Confirm: one Video, the attachment as first frame, 9:16, Lite, 720p,
   1 output, 8 generated seconds, action complete by second 6.
6. Read the displayed credit charge before clicking. It must describe one
   generation. If the interface forces two outputs, cancel and recheck the
   output setting instead of approving an accidental double generation.
7. Click **Generate**.
8. Download the result into `../renders/flow/drafts/` as
   `S02_flow_lite_720p_take02.mp4`.
9. Watch it at normal speed and 0.5× speed.
10. Complete the S02 row in `../REVIEW_SHEET.csv`.
11. If it passes, stop generating S02 and move to S05.
12. If it fails, choose one repair from `../PRODUCTION_WORKFLOW.md`, change only
    that sentence, and generate one more Lite output.
13. If the same defect remains, stop using Lite for this shot. Try one Fast
    output or one Kling Standard fallback.

Repeat with the next row. Remove the previous start frame and prompt before
loading the next shot.

## Prompt-first repair rule

Do not use synthetic zoom, freeze frames, warping or local motion effects to
turn a failed generated shot into accepted footage. A normal duration trim is
part of timeline assembly, but it may only select an already clean continuous
segment. If the camera, fluid motion or rigid cube fails, record the defect in
`../REVIEW_SHEET.csv`, revise one specific prompt instruction and generate a
new take when credits are available.

## When to use Fast

Use `Veo 3.1 Fast` only when:

- the Lite composition and endpoint are correct;
- a motion or detail defect remains visible at phone size;
- the full 720p rough cut shows that this particular shot needs improvement.

Keep 720p, one output, the same start frame, and the accepted prompt. Append
`_fast_take01.mp4` to the filename.

## When to use Quality

Quality is optional. Consider it only after the complete rough cut passes, and
only for S02 or S06 when an A/B test shows a visible phone-sized improvement.
Generate one output, trim it in the editor, and do not run Quality across all
six shots by default.

## Do not proceed when

- the cube deforms;
- more than one coral filament appears;
- unreadable AI writing or an equation appears;
- the camera misses the requested endpoint;
- fewer than 8–12 clean frames exist at either end.

Do not open Storyboard Studio to fix these problems. In the expanded Agent
fallback, remain in the same session and use only the exact targeted retry
block; do not ask the Agent to invent a storyboard or multiple alternatives.
