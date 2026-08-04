# Edit, motion and sound plan

## Sequence setup

- Canvas: 1080 × 1920, 30 fps, Rec.709, progressive.
- Source scouts: 720 × 1280. Scale them to fill the vertical canvas; do not
  upscale individual shots before the complete phone rough cut passes review.
- Duration: 00:29.70 maximum.
- Keep at least 120 px from left/right edges, 210 px from the top, and 350 px
  from the bottom for critical mid-film copy.
- S01, the first question in S02 and S07 are deliberate full-page exceptions.
- Generate footage without captions; use the SRT and vector overlays in post.

## Exact research-clock mapping

- Research window: `00:12.40–00:25.20` = 12.80 s.
- Real interval: `T+00:00:00–T+33:28:21`.
- Scale: about 0.382 seconds of film per research hour.
- Put the persistent `RESEARCH TIME T+hh:mm` clock near the upper safe edge.
- Animate the clock linearly. Snap event chips only at the film times in
  `TIMELINE_EVIDENCE.md`.

Do not ease or pause the clock. The feeling of slowdown comes from shot design:
S04 has many rapid cuts, S05 has fewer and longer inspections, and S06 uses one
slow move with a long terminal hold.

## Shot retiming

### S01 — 00:00.00–00:03.50

- 12-frame fade from black.
- 100% stable plate, 4% total dolly.
- `I ASKED AI TO ATTACK MATH` appears at 00:00.30; reveal `I BARELY
  UNDERSTAND` at 00:01.00.
- Add `HUMAN-DIRECTED • AI-EXECUTED` and `EVIDENCE, NOT TRUST` as quiet footer
  lines at 00:01.65.
- Keep the value proposition readable through 00:03.30.

### S02 — 00:03.50–00:09.50

- 03.50–04.20: 100%.
- 04.20–05.80: ramp to 450%.
- 05.80–06.80: ramp down to 35% on the coral bend.
- 06.80–09.50: 70% macro hold for the two labels and pointer.

### S03 — 00:09.50–00:12.40

- Use a restrained macro orbit, no speed ramp.
- Build one five-square loop: `PROOF GAP`, `ARTIFACTS`, `SEARCH`, `AI
  ANALYSIS`, `NEW LEMMA ↺`.
- Reveal the complete loop by 00:10.15, then pulse through the nodes once. Keep
  `METHOD KINSHIP • PROTAS / TAO` as a small qualifier, not a headline.
- Do not show portraits; the visual subject stays the fluid and the method.

### S04 — 00:12.40–00:17.21: fast early phase

- Use an 8-second shell-flight plate.
- Retiming curve: 260% → 800% → 320%.
- Each early event uses: 5-frame white candidate, 4-frame amber pulse,
  10-frame coral hold, a 4-frame teal refinement pulse, then exit.
- Keep `AI PROPOSES → C++ TESTS → SAVED EVIDENCE → AI REFINES ↺` as a thin
  process rail and `CHEAP GATES FIRST` as its readable headline.
- Land chips at 12.63, 14.89, 15.91 and 17.21.
- The F002 mark may be almost subliminal; F012 gets the longest hold.

### S05 — 00:17.21–00:21.04: deeper optimization

- Keep the base plate between 55% and 80% speed.
- Change the stage headline to `SURVIVED → DEEPER SEARCH`.
- F013: 14-frame hold at 19.89; F014: 16-frame hold at 20.34.
- At 21.04, make a one-frame black dip plus dry click for the memory failure.
- Resume on the identical composition with `CHECKPOINT → RESUME`; no
  glitch montage or fake terminal output.
- Use one thin pointer from the active card to the scan/filament intersection.

### S06 — 00:21.04–00:25.20: hardest phase

- Slow the generated pullback to 40–55%.
- Introduce `AI-REFINED FRONTIER` before the K12 context line.
- Introduce `K12 • 4.88 GiB` quietly; it is context, not the headline.
- Hold PNT-13 rejection for 18 frames at 24.86.
- At 25.03, arrest almost completely on `PNT-12 • SURVIVED THE TEST`.
- Add `THIS IS NOT A PROOF` immediately; no glow, checkmark or resolved chord.

### S07 — 00:25.20–00:29.70

- Use the wide final frame at 98–100% speed with only a 2% pullback.
- Place a 72% near-black veil over footage for contrast.
- Reveal the results as one group at 25.35 and the estimate group at 27.10.
- Lead with `AI-EXECUTED CAMPAIGN • 33 h 28 min` and close with `NOT A PROOF •
  EXPERT REVIEW REQUIRED`.
- Hide regular SRT subtitles in this shot because the card itself carries the
  final words.
- Keep `Clay solution: no honest estimate` visible for at least 1.8 seconds.

## Pointer line

- Card anchor: approximately `(x=180, y=420)`.
- Target: the amber scan / coral filament intersection in the chosen crop.
- Draw a 2 px off-white Bezier line; animate it over 6 frames.
- Finish with an 8 px hollow coral circle. It must point to the visible bend.

## Status overlays

| Status | Color | Treatment |
| --- | --- | --- |
| AI-proposed candidate | White | Quiet card, no pulse |
| Deterministic C++ test | Amber | One short scan pulse |
| Rejected | Coral | Dry stamp + short hold |
| Evidence saved / AI refines | Teal | Ledger mark + short forward pulse |
| Proved partial estimate | Teal | Stable underline, no glow burst |
| Survived finite test | White + amber edge | Always paired with `not proof` |
| Open | White + cyan edge | Open-ended line continuing right |

## Sound design

- Bed: restrained 90–105 BPM electronic pulse, no heroic trailer rise.
- Flow: low glassy water movement and soft broadband whoosh.
- S02: one Doppler whoosh, then near-silence at the macro arrest.
- S04: tight dry ticks; let their frequency imply fast iteration.
- S05: widen the gaps between ticks; use one dry dropout for the crash/resume.
- S06: one low scan tone and a long airy tail.
- S07: no triumphant cadence; voiceover stays 5–6 dB above the bed.

## Mix and export

- Dialogue peaks: about −3 dBFS.
- Integrated loudness: approximately −14 LUFS; true peak ≤ −1 dBTP.
- Finish the complete 720p-source rough cut first. Upscale only the locked
  timeline, and only if the delivery platform or phone review justifies it.
- Master: ProRes 422 HQ or DNxHR HQX, 1080 × 1920, 30 fps.
- Social: H.264 High Profile, 15–25 Mbps, AAC 48 kHz / 320 kbps.
- Export burned-caption and clean-captionless versions.
- Thumbnail: S01 at 01.60 s, title no smaller than 92 px.
