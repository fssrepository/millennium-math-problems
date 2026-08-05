# Navier–Stokes video — compact session handoff

Updated: 2026-08-05

This is the resume point for a fresh **video-only** session. Do not reconstruct
the work from chat history or preload the full mathematical campaign. The
detailed production log remains available for audit, but this file routes the
minimum context needed for the next action.

## Scope

- Work only inside this video directory unless the user explicitly expands the
  task.
- Do not edit repository-level `AGENTS.md`, research code or proof artifacts.
- Do not read the archived `.codex` session by default. The verified facts
  needed by the film are already captured in `sources/RESEARCH_NOTES.md` and
  the story documents.
- Preserve the current numbered folder structure; it was introduced to keep a
  fresh session from loading unrelated material.

## Current folder map

- `00_START_HERE.md` — human-facing production entry point.
- `01_story/` — story, scientific evidence, visual truth and narration.
- `02_production/` — workflow, timing, review decisions and durable log.
- `03_platforms/` — Flow, Higgsfield and OpenArt prompts/runbooks.
- `04_editorial/` — overlays, captions and release copy.
- `assets/`, `renders/`, `scripts/`, `sources/` — stable technical locations.

## Film and claim lock

- Format: 9:16, 29.70 seconds, English narration/captions.
- The film shows **human-directed, AI-executed** research, not autonomous proof.
- The Millennium problem is not solved; expert review is still required.
- Sixteen candidate routes were rejected for different reasons.
- Only F015/PNT-13 receives a detailed visual explanation: separated shell
  correlations were expected to decay but stayed near one. The evidence is
  saved and the active frontier changes to the still-open PNT-12 target.
- The visual is a concept visualization, not the C++ program's numerical output.

## Current production checkpoint

### Latest checkpoint — supersedes the pre-generation action below

- V03's widened dissolve still could not make two persistent geometries read
  as one. Exact frame matching identified approximately one-second overlap
  loops: source frame 359 recurs at 383, and source frame 551 recurs at 575.
  The current file is
  `renders/edit/navier-stokes-flow-post-v04-review-720p.mp4`. It removes those
  two actual repeated intervals and gently retimes each remaining continuation
  to the same eight-second slot. There is no visual dissolve at either join;
  boundary contact strips and full decode passed.

- V02's two-frame bridge at the 15-second boundary still appeared to jump in
  user review. The current file is
  `renders/edit/navier-stokes-flow-post-v03-review-720p.mp4`. It retains only
  two dropped source frames but expands that single bridge to six 30 fps output
  frames from 14.80 to 15.00. Frame-strip and full-decode QC passed. The
  23-second bridge remains the accepted two-frame version.

- User QC superseded the 1080p post v01: another geometry jump was visible
  around 13–14 seconds and its de-clicked audio sounded distorted in several
  places.
- The current review is
  `renders/edit/navier-stokes-flow-post-v02-review-720p.mp4`: native 720 ×
  1280, 30 fps, 891 frames and 29.700 seconds. It repairs the persistent source
  jumps at exactly 15.000 and 23.000 seconds with two dropped source frames
  plus a two-frame dissolve at each boundary.
- V02 removes broadband de-click entirely. Only a 12 ms seam fade at 7 seconds
  and 83 ms equal-power crossfades at 15 and 23 seconds touch the audio. The
  remainder is the original Flow score with a safe gain adjustment.
- Review V02 before any upscale. If an audible defect remains away from those
  joins, record its exact time because it originates in the Flow master and
  needs a separate local decision.

- The paid Flow chain is complete: V01 plus three Extend actions used 40
  credits. The held fifth 10-credit slot remains unspent.
- The downloaded continuous source is
  `/home/raxim/Downloads/Initial_Scene_-_2026-08-05_202608051411.mp4`: 720 ×
  1280, 24 fps, 31.000 seconds with generated stereo audio.
- The first caption-led post master is
  `renders/edit/navier-stokes-flow-post-v01-29.70s.mp4`: 1080 × 1920, 30 fps,
  891 frames and exactly 29.700 seconds. It decodes without error.
- The only meaningful visual join defect was the user-observed jump at source
  23 seconds. Post drops two source frames and applies only a two-frame
  dissolve. The large S06 frequency-shell card naturally reduces the remaining
  visibility without inventing a new scene.
- The generated score is preserved. De-click is local to the first two joins;
  the 23-second sound follows the same short crossfade as picture. Measured
  peak is −1.11 dBTP.
- Exact opening/topic copy, research overlays, separate tail diagram,
  expanding glass cards, border light sweep and final localized glints are now
  deterministic post effects. A spoken English voice has not been baked
  because no voice asset is selected; the master is caption-led.
- Next action: user phone review of that post master. Do not spend the reserve
  on a retry unless the defect survives local repair and materially harms the
  film.

The remaining bullets in this section document the earlier planning and
pre-generation checkpoint; they are retained for audit, not as the current
instruction.

The durable source of truth is `02_production/PRODUCTION_LOG.md`. Its present
checkpoint is:

- The user-selected three-second-per-frame reference slideshow exists at
  `renders/edit/navier-stokes-10-frame-morph-previs-v01-30s.mp4`. It was
  restored byte-for-byte from commit `46ac16a`; FFprobe verifies 720×1280,
  24 fps and exactly 30.000 seconds. Keep it for visual review, but do not use
  it to approve the still-missing PNT-13/PNT-12 story transition.

- A no-credit v02 story preview now exists at
  `renders/edit/navier-stokes-10-frame-story-previs-v02-30s.mp4`. FFprobe
  verifies 720×1280, 24 fps and exactly 30.000 seconds. K04–K09 now use
  progressively tighter crops of one source, so the far-shell diagram appears
  over the same forward-moving fluid plan rather than after a background-still
  swap. K06–K09 move deeper
  into the same fluid while editor overlays show expected decay → persistent
  far-shell response → saved PNT-13 rejection → PNT-12 open. The final three
  seconds now use a true scale pullback on the same source image: the same
  coral structure remains visible in the central cube while weather, flight,
  blood-flow and energy forms are revealed around it. It still requires user
  approval after the latest camera-path rebuild.

- The paid route is now one fresh eight-second Veo 3.1 Lite video followed by
  three sequential **Extend** actions on that same video. These are four fixed
  generation slots, not four scenes. Budget: 40 credits, with the fifth
  10-credit slot reserved for one targeted retry. The complete prompts are in
  `03_platforms/flow/CHAIN_32S_RUNBOOK.md`.

- Audio is now part of the same continuity contract. V01 generates one
  restrained premium science-documentary electronic score; V02–V04 explicitly
  continue its tempo, tonal center, instrumentation and ambience without a
  restart at the Extend boundaries. No separate music track is planned. Veo
  must not generate voice, speech, singing or lyrics; narration and captions
  remain editor-controlled.

- The static final preview strip illustrates only composition and camera
  scale. It is not the intended V04 motion. The paid V04 must keep simultaneous
  motion in the revealed field: rotating weather vortex, travelling wing
  streamlines, pulsing vessel flow, spinning turbine rotor and propagating cyan
  wavefronts. Reject frozen application scenery.

- S06 no longer points to any location in the liquid. The tail is shown only as
  the far end of the separate frequency-shell diagram, not as a physical wave
  tail or raw spatial field location. No arrow, target ring or connector into
  the fluid is allowed.

- Major captions use an editor-owned light animation: a compact dark-glass box
  enlarges, a restrained cyan/coral sweep crosses its border, and the essential
  lines appear quickly from top to bottom. In S07, short tracked glints accent
  the already-moving weather vortex, wing airflow, vessel flow and turbine.
  The exact paste-ready direction is `04_editorial/OVERLAY_MOTION_PROMPT.md` and
  must never be pasted into Flow/Veo.

- The final post pass identifies the subject immediately at the beginning with
  `NAVIER–STOKES REGULARITY` and `CAN A SMOOTH FLUID BECOME INFINITELY ROUGH?`.
  This is added over the clean V01 opening and does not require regeneration.

- The user reported refreshed credits. Create a new Flow project named
  `Navier–Stokes — Continuous Morph V01–V04 — 9x16`. Do not use the old
  `Lemma Stress Test` project's visible pending approval card.

1. Preflight review of the anchor, v02 story preview, prompt chain, audio
   contract, file paths and credit gate passed on 2026-08-05.
2. No generation credits should be spent ahead of the supervised steps.
3. Create the clean new Flow project and verify its settings before approving
   V01; then accept each result before moving to the next Extend.

Minimum files for that action:

- `01_story/FALSIFICATION_VISUAL.md`
- `assets/previs-10-v02/README.md`
- `03_platforms/flow/CHAIN_32S_RUNBOOK.md`
- `02_production/PRODUCTION_LOG.md` only when an exact prior decision or source
  filename is needed

Do not preload the rest of the bundle.

## Existing media decisions

- S01 raw Flow take exists and still needs full raw review.
- S02 take 02 is a usable provisional plate, not final-locked.
- S05 motion is a reference only because the generated amber beam/ring makes it
  unacceptable for the final film; local removal damaged the fluid.
- S06 take 01 is rejected: it rendered numeric labels and morphed markers.
- The still-animation S06 fallback is also rejected because it froze the fluid
  and the synthetic zoom looked unstable.
- Future S06 must have living fluid, a level single optical axis, and no numeric
  instructions or generated status graphics. Exact guides, labels, pointer and
  ledger state belong in post.

## Audit and cleanup rule

Keep `02_production/PRODUCTION_LOG.md` and `02_production/REVIEW_SHEET.csv`:
they record why a take was rejected and how the current route was reached.
Do not keep duplicate, rejected, failed-repair or superseded video binaries in
Git. A removed binary's filename and decision remain in the audit documents;
repository history must not be rewritten for cleanup.

## Working behavior for a fresh session

- Before a long operation or image generation, tell the user what is starting
  and report back regularly; do not disappear while a batch runs.
- Do not generate a large image pack in one call. Work checkpoint by checkpoint
  so a local revision normally affects only the adjacent material.
- Do not rebuild an animatic or spend credits merely to verify this handoff.
- Preserve unrelated worktree changes and never flatten the new directory map.

## Completion condition for the next checkpoint

The next checkpoint is complete after the user watches
`renders/edit/navier-stokes-flow-post-v04-review-720p.mp4` at phone size and either
accepts it or identifies an exact timecoded defect. Use deterministic post for
any local typography, mix or join correction. Spend the held 10 credits only
if a material generative defect cannot be repaired locally and the user
explicitly approves that targeted retry.
