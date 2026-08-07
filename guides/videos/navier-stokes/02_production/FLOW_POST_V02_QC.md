# Flow post v02 review — quality-control record

Reviewed: 2026-08-05

- File: `../renders/edit/navier-stokes-flow-post-v02-review-720p.mp4`
- SHA-256: `a76f9e7a224dabefae29845467fabeb49c7608844d57cd517dfc13c4bd978acc`
- Video: H.264, 720 × 1280, 30 fps, 891 frames
- Runtime: 29.700 seconds
- Audio: AAC stereo, 48 kHz; approximately −16.28 LUFS and −1.33 dBTP
- Full decode: passed

## Changes from v01

- Kept review resolution at the 720 × 1280 source size; no upscale.
- Added the same minimal repair at the 15.000-second source jump as at the
  23.000-second source jump: drop two frames, dissolve for two frames.
- Removed full-band de-click processing.
- Limited audio work to a 12 ms fade at the 7-second join and 83 ms
  equal-power crossfades at the 15- and 23-second joins.

## Human gate

Listen through once and report any remaining audible defect with an exact
timecode. A defect away from 7, 15 or 23 seconds is part of the generated Flow
score and must not be treated with an indiscriminate whole-track filter.
