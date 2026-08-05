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

The durable source of truth is `02_production/PRODUCTION_LOG.md`. Its present
checkpoint is:

- The user-selected three-second-per-frame reference slideshow exists at
  `renders/edit/navier-stokes-10-frame-morph-previs-v01-30s.mp4`. It was
  restored byte-for-byte from commit `46ac16a`; FFprobe verifies 720×1280,
  24 fps and exactly 30.000 seconds. Keep it for visual review, but do not use
  it to approve the still-missing PNT-13/PNT-12 story transition.

1. No additional video-generation credits should be spent yet.
2. The superseded ten-frame v01 morph did not visibly explain the PNT-13
   falsification and must not approve the story.
3. Rebuild/composite planning states K06–K09 so the v02 previs visibly shows:
   expected decay → persistent far-shell response → saved PNT-13 rejection →
   PNT-12 open frontier.
4. Review that v02 at phone size against every acceptance item in
   `01_story/FALSIFICATION_VISUAL.md`.
5. Only after user approval may the Flow/OpenArt/Higgsfield prompts be revised
   and generation resume from the second video piece.

Minimum files for that action:

- `01_story/FALSIFICATION_VISUAL.md`
- `assets/previs-10/README.md`
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

The next checkpoint is complete only when a phone-size v02 planning preview
clearly communicates all four PNT-13/PNT-12 states and the user approves it.
Prompt rewrites and new paid generations remain downstream work.
