# Keyframe assets and provenance

These project-bound reference images were generated with Codex's built-in
`imagegen` workflow on 2026-08-04. They are concept visuals, not outputs of the
Navier–Stokes C++ laboratory.

| File | Role | Generation method |
| --- | --- | --- |
| `keyframe-01-world.png` | Opening / master world identity and corrected S02 start frame | New 9:16 generation, then precision image edit after the first Flow scout made the coral feature behave like a solid braided rope |
| `keyframe-03-checkpoint.png` | Interior macro / falsification focus | Derived from the master; cable-like coral material replaced with restrained fluorescent dye |
| `keyframe-06-frontier.png` | Closing shell/ledger world | Derived from the master; corrected to 16 marks, then central coral material corrected without changing the ledger count |
| `keyframe-06-frontier-heightfield-v2.png` | Rejected nested-cube replacement concept | Generated from the canonical frontier; uses curved cyan height-field ridges but still contains a baked amber path |
| `keyframe-06-frontier-heightfield-clean-v3.png` | Clean frontier previs base | Precision edit of v2 with the baked amber path removed for deterministic progress animation |
| `overlay-top-gradient-1080x1920-rgba.png` | Source overlay for darkening bright particles behind upper-screen copy | Deterministic transparent SVG gradient rendered as a true RGBA PNG; no generative model |

## Continuous Flow chain anchors

`flow-chain/` contains exact endpoint frames for the three-generation,
32-second continuity route documented in `../flow/CHAIN_32S_RUNBOOK.md`.

| File | Provenance and role |
| --- | --- |
| `flow-chain/chain-anchor-00-s01-last.png` | Exact final frame of accepted raw S01; first frame of C01 |
| `flow-chain/chain-anchor-01-s02-last.png` | Exact final frame of raw S02 take 02; last frame of C01 and first frame of C02 |
| `flow-chain/chain-anchor-02-s05-last-clean.png` | Exact clean final frame of raw S05; last frame of C02 and first frame of C03 |
| `flow-chain/chain-anchors-contact-sheet.png` | Left-to-right check of the three extracted anchors and canonical `keyframe-06-frontier.png` endpoint |

The original generated clips remain preserved for audit. These PNGs are
unmodified extracted source frames, not locally reconstructed animation.

`previs-10/` contains the current ten-frame, three-second-step visual plan and
its own provenance table. It is the source for the 30-second morph preview.

The reusable text prompt is in `../CREATIVE_BIBLE.md` under **Master visual
prompt** and **Shared negative prompt**. The more detailed motion prompts are in
the three platform `PROMPTS.md` files.

Visual invariants:

- one glass cube;
- one teal fluid system;
- one recognizable coral vortex filament;
- amber only for testing/checkpoint accents;
- no generated typography;
- no “problem solved” imagery.

The canonical `keyframe-01-world.png` now contains the corrected compact
fluorescent-dye vortex: it has diffuse internal flow wisps, no solid braided
surface, and visible clearance from every tank boundary. The rejected first
Flow scout records how the earlier source behaved; do not re-upload a locally
cached copy of the earlier image.

The same material audit was applied to the canonical S03/S05 macro and S06
frontier frames. `keyframe-03-checkpoint.png` still contains the earlier amber
scan intersection and is now a legacy composition reference only; the final
S05/S06 route must contain no generated scan line or beam. The deterministic
expected-decay guide, measured far-shell response, pointer, saved rejection and
PNT-12 marker are defined in `../FALSIFICATION_VISUAL.md` and added in post.
`keyframe-06-frontier.png` was visually recounted and still contains exactly 16
separate checkpoint cubes.

`overlay-top-gradient-1080x1920-rgba.png` is a compositing source, normally
placed above the generated plate and below typography. Its alpha is already
embedded. CapCut Online flattened both the indexed-alpha and true-RGBA tests to
an opaque black card on 2026-08-04, so do not use this as a separate CapCut
track. Use a native editor gradient or the precomposited S02 file documented in
`../flow/RUNBOOK.md`. Desktop compositors that preserve PNG alpha may leave the
overlay opacity at 100%. The editable source is the adjacent SVG.
