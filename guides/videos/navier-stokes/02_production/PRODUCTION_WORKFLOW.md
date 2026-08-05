# Cost-efficient manual video-AI production workflow

This workflow spends only after the previous result has passed review. It uses
720p scouts because the film is designed for a phone screen. A higher-quality
generation is optional and justified shot by shot, never assumed for all six.

## Spending logic

```text
Free still-image animatic
  ↓ approve timing and copy
One 720p Lite scout for S02
  ↓ inspect
One 720p Lite scout for S05
  ↓ inspect
One 720p Lite scout for S06
  ↓ approve visual world and camera grammar
One 720p Lite scout each for S01, S03, S04
  ↓ assemble complete mobile rough cut
Repair only failed shots
  ↓ phone review
Optional Fast or Quality upgrade for visible defects only
```

Do not generate multiple variants by default. One output is a test; the next
generation is authorized only by a concrete finding from that test.

## Model ladder

1. **Flow / Veo 3.1 Lite, 720p** — default scout for every shot when free
   account access is available.
2. **Flow / Veo 3.1 Fast, 720p** — only when the Lite composition works but
   motion or fidelity is visibly insufficient.
3. **Higgsfield or OpenArt / Kling 3.0 Standard, 720p** — only when Flow cannot
   execute a required camera move after one targeted retry.
4. **Veo 3.1 Quality** — optional for at most one or two hero shots after the
   complete rough cut passes on a phone. Skip it when the 720p take already reads.

Use only one paid platform at a time. Do not buy OpenArt and Higgsfield access
in advance; choose one fallback only after a specific Flow failure identifies
what is missing.

Google's current help pages both describe 50 daily credits for eligible
non-subscribers and, elsewhere, subscription-gated Flow access. Test the actual
account before treating Flow as free. If access is blocked, skip directly to
one paid Kling route; do not buy Google access merely to preserve this ladder.

## No plugin is required

The full film can be made with Flow's standard prompt box, Video → Frames, and
a conventional editor. Do not use Agent, Storyboard Studio, Smart Shot, AI
Director, Motion Control, Cast, Soul ID, or a third-party plugin.

## 0. Zero-generation animatic

Before opening a video generator:

1. Create a 1080 × 1920, 30 fps editor timeline.
2. Place the ten planning frames from `../assets/previs-10/` at their three-second
   story intervals.
3. Set shot lengths from `shot-manifest.csv`.
4. Add the copy from `../04_editorial/EDITOR_OVERLAYS_EN.md`.
5. Import `../04_editorial/captions/navier-stokes-en.srt` and a rough narration from
   `../01_story/narration/en.md`.
6. Watch the still-image cut on a phone.

Do not spend generation credits until the title, narration, final estimate,
safe margins, and 29.70-second pacing work in this animatic.

## 1. Generate the three high-risk scouts first

Generate in this order:

1. **S02** — fastest push-in and macro endpoint;
2. **S05** — stable deeper-optimization motion with no scan line or beam;
3. **S06** — continuous fluid inspection supporting the PNT-13 obstruction
   and PNT-12 handoff defined in `../01_story/FALSIFICATION_VISUAL.md`.

For each shot:

1. Select 720p, 9:16, Veo 3.1 Lite, and one output.
2. Upload the assigned start frame.
3. Paste the shot's complete prompt.
4. Check that the interface's displayed charge is for one generation. If it
   shows multiple outputs, cancel and correct the output setting.
5. Generate once and download immediately.
6. Watch at normal speed and 0.5× speed.
7. Fill in `REVIEW_SHEET.csv`.
8. If it passes, stop generating that shot and continue to the next one.
9. If it fails, identify one defect and make one targeted Lite retry.
10. If the same defect remains, stop spending on Lite and use one Fast attempt
   or one Kling 3.0 Standard fallback—not both simultaneously.

Do not generate S01, S03, or S04 until these three shots establish that the
cube, filament, liquid material, and camera language are viable.

## 2. Complete the first 720p pass

Generate one Lite scout each for:

1. S01;
2. S03;
3. S04.

Use the same review gate. Do not generate an alternate take when the first take
passes. S07 costs nothing: it is a hold on the fixed S06 endpoint image.

## 3. Mandatory review gate

Accept a shot only when every relevant condition passes:

- the cube stays rigid and neither melts nor breaks;
- exactly one coral vortex filament remains;
- the teal flow looks like liquid, not smoke or fabric;
- no generated text, watermark, person, or pseudo-equation appears;
- the camera stops at the requested final composition;
- at least 8–12 clean edit frames exist at both ends;
- the shot-specific requirement below passes.

| Shot | Most important acceptance condition |
| --- | --- |
| S01 | The upper 38% is calm and dark enough for the title |
| S02 | The glass stays intact and the coral bend stops at center-right |
| S03 | The macro bend remains in focus; no person or portrait appears |
| S04 | Four distinct test pulses and stable shell geometry remain visible |
| S05 | The coral filament stays stable and no vertical scan line or beam appears |
| S06 | One continuous fluid state supports expected decay, a persistent far-shell response and the PNT-12 handoff |

## 4. One-defect retry table

| Problem | Single change for the next attempt |
| --- | --- |
| Coral filament duplicates | Add: `one single continuous coral filament, never duplicate` |
| Cube deforms | Add: `rigid unchanged glass cube`; reduce camera speed |
| Flow looks like smoke | Add: `coherent heavy liquid, not smoke, not vapor` |
| Camera misses target | Keep model and image; revise only the camera sentence |
| Motion is excessive | Reduce camera speed; if it repeats, use Kling Standard |
| Generated writing appears | Start the negative block with: `absolutely no visible text or symbols` |
| S06 invents markers, percentages or diagrams | Generate clean fluid motion only; build the separate decay card, ledger and PNT-12 marker in post from `../01_story/FALSIFICATION_VISUAL.md`; do not connect the card to the fluid with a pointer |

Never change the model, start frame, and prompt together. That makes the failed
generation diagnostically useless.

### Logged S02 source-repair exception

Take 01 exposed two independent causes that were visible in the evidence: the
master still already described the coral feature as braided, while the phrase
`accelerates through the intact front glass` encouraged a physical impact. For
take 02, keep Veo 3.1 Lite, 720p, duration and output count fixed, but use the
corrected canonical start frame and the non-contact optical-camera
sentence in `../03_platforms/flow/AGENT_READY_TO_PASTE.md`. This compound repair is intentional
and logged; it avoids spending another 10 credits on either already-confirmed
defect. Resume the one-defect retry rule after take 02.

## 5. Assemble before upgrading

1. Place the six accepted 720p scouts in a 1080 × 1920 editor timeline.
2. Cut the sequence to exactly 29.70 seconds.
3. Apply the speed and freeze points in `EDIT_PLAN.md`.
4. Add the English overlays, SRT, rough narration, and temporary sound bed.
5. Composite the S06 expected-decay guide, measured far-shell response,
   saved-rejection ledger and orange PNT-12 marker exactly as specified in
   `../01_story/FALSIFICATION_VISUAL.md`. Keep the diagram separate from the
   liquid with no pointer or target ring.
6. Watch the full rough cut on the target phone at normal brightness.

Only a defect visible in this complete phone review justifies another video
generation. A flaw visible only at 200% desktop zoom does not.

## 6. Selective finishing

- If every 720p scout reads cleanly on a phone, keep them.
- If one shot lacks motion quality but has the right composition, regenerate
  only that shot with Veo 3.1 Fast.
- If a camera path is still wrong, use one Kling 3.0 Standard attempt on one
  fallback platform.
- Consider Quality only for S02 or S06, and only when an A/B comparison shows a
  visible improvement at phone size.
- Upscale or export to 1080 × 1920 only after picture lock. Do not upscale
  rejected drafts.

## 7. Final acceptance

Complete `DELIVERY_QC.md`, then watch the export:

1. on a phone with sound;
2. on a phone while muted;
3. at 0.75× speed;
4. once while reading only the captions.

Export only after all four reviews pass.
