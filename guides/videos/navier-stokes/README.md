# Navier–Stokes short-video production bundle

Status: production guide for a conceptual 9:16 explainer, not a numerical simulation export  
Generation target: 720 × 1280 first; optional upscale only after picture lock
Timeline target: 1080 × 1920, 30 fps, 29.7 seconds, English narration and captions

## What story this bundle tells

The film is not simply a fluid simulation or a list of numerical results. It
shows a human-directed, AI-executed research workflow. The human supplied the
question, constraints, laptop and continuation decisions. Codex wrote and ran
the C++, inspected failures, designed narrower lemma candidates and stored the
evidence needed to continue after a crash. The result is a reproducible proof
frontier for expert review—not an automated proof and not a Clay solution.

The five-square loop visible in S03 is the governing structure:

`PROOF GAP → ARTIFACTS → SEARCH → AI ANALYSIS → NEW LEMMA ↺`

The reason for building the repository is practical: a cheap, deterministic
obstruction can eliminate a weak direction before a researcher spends weeks or
months trying to prove it. More compute is reserved for the candidates that
survive the early gates.

## Recommended route

For this exact film, the strongest workflow is:

1. Build the complete 29.7-second still-image animatic before spending credits.
2. If the account exposes free Flow generation, test only S02, S05, and S06 at
   720p with one **Veo 3.1 Lite** output each. Otherwise choose one paid Kling
   fallback from `02_production/PLATFORM_COST_GUIDE.md`.
3. If those pass, generate S01, S03, and S04 one at a time with Lite.
4. Use **Veo 3.1 Fast** only for a shot that fails visually at Lite quality.
5. Use **Veo 3.1 Quality** only for one or two hero shots if a phone review
   shows a visible need; it is not a mandatory final pass.
6. If a camera still fails after one targeted retry, remake only that shot with
   **Higgsfield or OpenArt basic image-to-video + Kling 3.0 Standard at 720p**.
7. Add every title, caption, arrow, freeze, time-ramp and counter with the
   reproducible `scripts/` pipeline or a normal editor such as DaVinci Resolve,
   Premiere or CapCut. Do not ask a video model to draw the English text.

Why this split: the film has no human character, so Flow/Veo's atmosphere and
scene fidelity matter more than face performance. Higgsfield/Kling is the
retake path for camera moves. OpenArt is an optional basic image-to-video
alternative. No Agent, Storyboard Studio, Smart Shot, Director, Motion Control,
character system or third-party plugin is required.

## Model choice at a glance

| Platform | Cheapest useful test | Selective upgrade | Best use here |
| --- | --- | --- | --- |
| Google Flow | Veo 3.1 Lite, 720p, one output | Fast only when needed; Quality only for a visible hero-shot gain | Main fluid footage |
| Higgsfield | Kling 3.0 Standard image-to-video, 720p, one output | Higher mode only after a successful scout | Camera rescue for S02/S05/S06 |
| OpenArt | Kling 3.0 Standard image-to-video, 720p, one output | Higher quality only after a successful scout | Complete fallback or alternate composition |

The recommendation is dated 2026-08-04. Model menus change, so the platform
guides include a fallback rule: preserve the supplied reference frames and
choose the newest model that supports 9:16 plus image/start-frame guidance.

## Folder map

- `00_START_HERE.md` — the first file to open; a plain-language folder map
  and seven-line workflow.
- `01_story/` — creative bible, storyboard, scientific evidence and narration.
- `02_production/` — workflow, edit plan, cost guide, manifests, review sheet
  and production log.
- `03_platforms/` — Flow, Higgsfield and OpenArt runbooks and prompts. Each
  platform has literal field-by-field instructions.
- `04_editorial/` — English overlays, captions, manual/social copy and timed
  ASS overlay files.
- `assets/` — stable reference frames, overlays and previs planning images.
- `renders/` — prepared local drop zones for drafts, finals, audio and edits.
- `scripts/` — deterministic FFmpeg repairs and the complete animatic build.
- `sources/RESEARCH_NOTES.md` — repository/session evidence and current product
  sources.

Each numbered directory contains its own short `README.md` index. Start with
`00_START_HERE.md`; use this README for the scientific and creative context.

## Fastest way to make it

1. Open `00_START_HERE.md`.
2. Run `scripts/build_animatic.sh` and watch the complete 29.70-second result at
   `renders/edit/navier-stokes-animatic-v02-raw-plates-720x1280.mp4` on a phone.
3. Mark each provisional plate or still placeholder `keep`, `repair` or
   `regenerate` before spending another credit.
4. Replace only the selected placeholders according to
   `01_story/STORYBOARD.md` and `02_production/EDIT_PLAN.md`.
5. Import the SRT, record or synthesize `01_story/narration/en.md`, then run
   `02_production/DELIVERY_QC.md`.

## Non-negotiable scientific wording

- The process is **human-directed and AI-executed**. Do not describe it as a
  fixed simulation, and do not describe it as autonomous proof.
- Codex wrote and ran the C++, interpreted failures, refined candidates and
  saved artifacts; the human set the objective and directed continuation.
- Expert mathematical review is still required.
- The Millennium problem is **not solved**.
- Navier–Stokes equations already support practical fluid models, including
  weather-related flows, aircraft aerodynamics and blood flow. A regularity
  proof would settle whether the idealized 3D mathematical model can develop a
  singularity; do not claim it would instantly solve those application areas.
- The footage is a **concept visualization**, not the output of the C++ run.
- `F000`–`F015` means **16 rejected lemma candidates**.
- Those candidates failed for different reasons. The film explains F015 /
  PNT-13 as one representative obstruction: distant frequency-shell
  correlations stayed near one instead of decaying. Do not apply that visual
  explanation to F000–F014.
- A finite search can falsify a candidate or let it survive the tested set; it
  cannot prove the cutoff-independent L4 statement.
- The demonstrated result is a fast, reproducible lemma-falsification loop,
  one proved partial far-tail estimate, and a still-open L4 closure.
- The `~1–3 days` estimate applies only to the next PNT-12 decision, not to a
  proof and not to the Clay problem.
- That estimate is observed laptop pace. `Up to ~10×` is an engineering ceiling
  for independent sweeps on ten comparable cluster nodes, not an end-to-end
  benchmark or proof acceleration factor.
- Six of the original seven Clay problems remain unresolved, with US$1 million
  allocated to each. Prize consideration requires a qualifying publication,
  at least two years, general mathematical acceptance and then CMI review; an
  unrefuted forum post is insufficient.
- Protas/Tao is a methodological comparison, not provenance of a combined
  proof; use the exact wording in `01_story/STORYBOARD.md`.

## The short version of the story

`Outside story → fluid problem and why it matters → Clay context → AI research
loop → fast early failures → one visible PNT-13 obstruction → saved evidence →
PNT-12 open frontier → laptop/cluster estimate → official prize path.`

From 00:12.40 to 00:25.20 the on-screen research clock is a linear compression
of 33 h 28 min. Early failures arrive quickly; later frames breathe longer as
the searches move from scaling and triads to exact gradients and K12.
