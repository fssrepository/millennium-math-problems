# Reproducible local video scripts

Run these from any working directory; each script resolves the guide directory
relative to itself.

- `build_10_frame_previs.sh` builds the current primary planning artifact,
  `../renders/edit/navier-stokes-10-frame-morph-previs-v01-30s.mp4`, from ten
  keyframes spaced at three-second story intervals. It spends no Flow credits
  and remains a superseded continuity reference.
- `build_10_frame_previs_v02.sh` builds the story-corrected no-credit review
  preview, `../renders/edit/navier-stokes-10-frame-story-previs-v02-30s.mp4`.
  Its K06-K09 boards follow one camera move into the same fluid while editor
  overlays show expected decay, the persistent far-shell response, the saved
  PNT-13 rejection and the PNT-12 open state. K10 plans the continuous
  application morph.
- `build_anchor_morph_previs.sh` builds the earlier five-anchor, 11-second
  direction test. It is superseded by the ten-frame previs.
- `build_flow_chain_32s.sh` is a legacy First/Last-frame assembly helper. The
  current Flow route uses the final V04 extended video directly.

- `build_animatic.sh` builds the complete 720 × 1280, 30 fps, 29.70-second
  review video at
  `../renders/edit/navier-stokes-animatic-v02-raw-plates-720x1280.mp4`. It uses
  the complete raw S02 source across S02 plus the first two seconds of S03 and
  a normal-speed raw S05 segment. Canonical stills fill missing shots, and
  `../04_editorial/overlays/animatic-overlays.ass` supplies timed typography.
- `render_s02_accept.sh` reproduces the provisional six-second S02 plate.
- `render_s05_salvage.sh` reproduces the provisional 3.83-second S05 plate.
- `render_s06_frontier.sh` reproduces a rejected local S06 experiment for
  audit only. It must not be used as final footage because it freezes the fluid
  and its synthetic zoom reads as camera instability.
- `render_s06_remove_label.sh` preserves the complete raw S06 duration, motion
  and source audio while covering only the generated `4%/6%` label.
- `build_raw_flow_review.sh` assembles all complete eight-second Flow clips,
  their original audio and labelled still placeholders into
  `../renders/edit/navier-stokes-raw-flow-review-v01-44p2s.mp4`. It does not
  retime or freeze generated footage.
- `build_four_flow_clips.sh` is the current minimal review build. It concatenates
  only S01, S02, S05 and the label-cleaned S06, in that order, for exactly 32
  seconds. It retains source audio and adds no labels, placeholders, trim,
  retiming, freeze or transitions.

The animatic carries a silent AAC track by design. Narration, music and final
sound effects are added only after the timing and shot replacements are locked.
