# Render drop zones

Download generated files into these folders using the naming convention in
`../PRODUCTION_WORKFLOW.md`.

- `flow/drafts/` — one-at-a-time 720p Veo Lite scouts and targeted Fast retries;
- `flow/finals/` — optional selected upgrades only after the rough cut;
- `openart/` — OpenArt alternatives;
- `higgsfield/` — controlled camera retakes;
- `audio/` — narration, music and sound effects;
- `edit/` — preview and final exports.

The current reproducible review export is
`edit/navier-stokes-animatic-v02-raw-plates-720x1280.mp4`. It is intentionally committed
with this production pass: 891 frames at 30 fps, exactly 29.70 seconds, with a
silent AAC track and burned-in rough-cut labels. Rebuild it with
`../scripts/build_animatic.sh`.

`edit/navier-stokes-animatic-v01-720x1280.mp4` is a superseded audit comparison
that used locally repaired S02 and S05 derivatives. V02 is the review source
because it preserves raw generated motion at normal speed.

`edit/navier-stokes-raw-flow-review-v01-44p2s.mp4` is the current untrimmed
motion review. It contains the complete eight-second S01, S02, S05 and S06
generations with source audio, plus labelled still placeholders for missing
shots. S06 has only its generated percentage label removed.

The newer primary review is `edit/navier-stokes-four-flow-clips-v01-32s.mp4`.
It contains only the four original eight-second Flow sources in story order,
with source audio and no added labels, placeholders, trim, retiming or freeze.
S06 has only its percentage label removed. S05 retains its generated amber
beam/ring because local removal damaged the fluid; it is queued for a no-scan
regeneration instead.

The `.gitkeep` files preserve the empty directory structure. Large generated
MP4/WAV outputs should only be committed if that is an intentional repository
policy decision.
