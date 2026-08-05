# Production control

Use this directory to plan, generate, review and deliver the film:

1. `PLATFORM_COST_GUIDE.md` — choose the model and spending route.
2. `PRODUCTION_WORKFLOW.md` — generate, inspect, retry and accept.
3. `COPY_PASTE_MAP.md` — shortest field-by-field production checklist.
4. `shot-manifest.csv` and `EDIT_PLAN.md` — timing and edit assembly.
5. `REVIEW_SHEET.csv` and `PRODUCTION_LOG.md` — decisions and current state.
6. `DELIVERY_QC.md` — final verification.

Paths stored in `shot-manifest.csv` are relative to the video bundle root, not
to this directory. This keeps them compatible with the scripts in `../scripts/`.
The CSV preserves the legacy editorial shot map. Current paid generation is
one V01 clip plus three Extends; use `../03_platforms/flow/CHAIN_32S_RUNBOOK.md`.
