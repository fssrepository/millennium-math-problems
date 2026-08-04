# Manual video-AI production workflow

This is not a one-shot generation. Create six separate clean clips, inspect
each one, proceed only after it passes review, then build the pacing and English
information layer in a conventional editor.

## Recommended route

```text
Keyframes
  ↓
Flow / Veo 3.1 Fast: 4 drafts per shot
  ↓ visual review
Flow / Veo 3.1 Quality: 2 final candidates
  ↓ only if the camera fails
Higgsfield / Kling 3.0: regenerate the problem shot
  ↓
Edit + speed ramps + freeze frames
  ↓
English overlays + narration + SRT
  ↓
Final phone review
```

OpenArt is a complete alternative route and a fast composition laboratory. You
do not need to generate every shot on all three services.

## No plugin is required

The full film can be made with Flow's standard prompt box, Video → Frames, and
a conventional editor. Do not use Agent, Storyboard Studio, Smart Shot, AI
Director, Motion Control, Cast, Soul ID, or a third-party plugin. Platform
runbooks mention these features only to tell you not to select them.

## 0. Preparation

Open these together:

- this file;
- the selected platform's `RUNBOOK.md`;
- its `PROMPTS.md` or `READY_TO_PASTE.md`;
- the `assets/` directory;
- `REVIEW_SHEET.csv`.

Download every AI output immediately. Suggested names:

```text
S01_flow_fast_take01.mp4
S01_flow_quality_take01.mp4
S02_higgsfield_kling_take03.mp4
```

Store them in the matching `renders/` subdirectory. A project preview on the
service is not a substitute for a local download.

## 1. Shot cycle — repeat for S01–S06

### 1.1 Input

1. Check the shot's assigned start frame in the platform runbook.
2. Set the aspect ratio to 9:16 and use the listed duration.
3. Paste the shot's complete prompt into the named field.
4. If there is a separate negative-prompt field, put the negative block there;
   otherwise leave it at the end of the main prompt.
5. Never ask the video generator for captions, equations, or interface text.

### 1.2 Draft generation

1. Produce three or four variants from identical inputs.
2. Do not change the model, prompt, and image at the same time; otherwise you
   cannot identify which change helped.
3. Download every variant and watch each at normal speed and 0.5× speed.

### 1.3 Mandatory review gate

Accept the shot only when every relevant condition passes:

- the cube stays rigid and neither melts nor breaks;
- exactly one coral vortex filament remains;
- the teal flow looks like liquid, not smoke or fabric;
- no generated text, watermark, person, or pseudo-equation appears;
- the camera stops at the requested final composition;
- at least 8–12 clean edit frames exist at both ends;
- the object and crop required by the next shot remain available.

Record the result in a working copy of `REVIEW_SHEET.csv`. A failed shot does
not proceed to the Quality pass.

### 1.4 Fix one problem at a time

| Problem | Next action |
| --- | --- |
| The coral filament duplicates | Keep the start frame and add: `one single continuous coral filament, never duplicate` |
| The cube deforms | Add: `rigid unchanged glass cube`; reduce camera speed |
| The flow looks like smoke | Add: `coherent heavy liquid, not smoke, not vapor` |
| The camera misses its target | Keep model and image; revise only the camera sentence |
| Motion is excessive | Request a slower camera or use the Higgsfield fallback |
| Generated writing appears | Start the negative block with: `absolutely no visible text or symbols` |
| S06 changes the 16 marks | Use generated motion only for the pullback, then cut to the fixed `keyframe-06-frontier.png` |

Make no more than two targeted retries on the same model. If the camera remains
unstable on the third result, move that shot to Higgsfield/Kling.

### 1.5 Final generation

1. Preserve the accepted start frame, prompt, and seed when available.
2. Switch to the platform's higher-quality mode.
3. Produce two final candidates.
4. Repeat the full review gate.
5. Choose the winning file before proceeding to the next shot.

## 2. Shot-specific gate

| Shot | Most important acceptance condition |
| --- | --- |
| S01 | The upper 38% is calm and dark enough for the title |
| S02 | The camera crosses the image plane without breaking the glass and stops on the coral bend at center-right |
| S03 | The macro bend remains in focus; no person or portrait appears |
| S04 | This is the fastest camera move; four distinct test pulses and stable shell geometry remain visible |
| S05 | This is slower and more precise; the amber scan actually crosses the coral filament |
| S06 | This is the slowest; an open teal path and a long stable final frame remain |

## 3. Editing workflow

1. Place the accepted clips in the order listed by `shot-manifest.csv`.
2. Cut the sequence to exactly 29.70 seconds.
3. Apply the speed and freeze points in `EDIT_PLAN.md`.
4. Counterexamples, `REJECTED` stamps, pointers, and the T+ clock are editor
   graphics; never regenerate them inside the clip.
5. At the end of S06, cut or dissolve to `keyframe-06-frontier.png` so the 16
   marks remain exact.
6. S07 is a 4.5-second hold on that fixed image with a 72% dark veil.

## 4. Text and audio

1. Copy the on-screen text from `EDITOR_OVERLAYS_EN.md`.
2. Import `captions/navier-stokes-en.srt`.
3. Record or synthesize the script in `narration/en.md`.
4. Target a narration duration of 29.2–29.7 seconds.
5. Hide ordinary subtitles during S07 because the full-page card already
   contains the result and estimate.

## 5. Final acceptance

Complete `DELIVERY_QC.md`, then watch the export:

1. on a phone with sound;
2. on a phone while muted;
3. at 0.75× speed;
4. once while reading only the captions.

Export the master and social version only after all four reviews pass.
