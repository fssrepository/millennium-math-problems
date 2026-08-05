# Render drop zones

Download generated files into these folders using the naming convention in
`../PRODUCTION_WORKFLOW.md`.

- `flow/drafts/` — one-at-a-time 720p Veo Lite scouts and targeted Fast retries;
- `flow/finals/` — optional selected upgrades only after the rough cut;
- `openart/` — OpenArt alternatives;
- `higgsfield/` — controlled camera retakes;
- `audio/` — narration, music and sound effects;
- `edit/` — preview and final exports.

Git stores only active selected source clips, accepted derivatives and deliberate
final deliverables. Rejected takes, local repair failures, review compilations,
morph previs, superseded animatics and backup exports are not retained as
binaries. Their decisions and defects remain auditable in `../REVIEW_SHEET.csv`
and `../PRODUCTION_LOG.md`.

Current retained motion assets:

- `flow/drafts/S01_flow_lite_720p_take01.mp4` — active opening source;
- `flow/drafts/S02_flow_lite_720p_take02.mp4` — selected provisional source;
- `flow/edits/S02_flow_lite_take02_6s_gradient_freeze.mp4` — active edit-safe
  derivative;
- `flow/drafts/S05_flow_lite_720p_take01.mp4` — provisional motion reference
  until its no-scan replacement exists.

Review and animatic scripts may write local outputs into `edit/`; those exports
must not be committed merely as backups. Commit a video only when it is the
current intentional review deliverable or final master. The `.gitkeep` files
preserve the empty directory structure.
