# Restorable Project history

Complete Project; 19 tasks. [Live Project](https://github.com/users/fssrepository/projects/5) · [Canonical complete backup](https://github.com/fssrepository/myscoutee-roadmap/tree/master/guides/statistics/math).

`project.json` preserves full issue bodies, statuses, typed field values and views. `measurements.json` preserves model/token/cost/time components and assumptions. `commit-map.json` links original and task-prefixed commit IDs. CSV/JSON provide convenient presentation exports. Recorded root intervals are retained in compressed JSON for timing reconstruction; there are no private chat texts.

## Latest measured task

MATH-19 is Done: **14,834,641 tokens**, **0.823521 active / 0.798249 allocated h**,
**$21.07 estimated Standard API equivalent**, **gpt-6-astra / high**, 141 responses
through `2026-10-25T16:43:58.451Z`. The complete session includes initial PNT-2
work and the continuation; the earlier 18 minutes covered only the continuation.

Current Project allocation: **1,300,786,118 tokens**,
**44.253883 allocated h**, **$748.75 API equivalent**.
The historical reconstruction is preserved separately. Account and task totals
overlap and are never added.

Reproduce usage from retained local history in the roadmap checkout:
`python3 guides/statistics/math/measure-math19.py`.

## Validate offline

Run from the repository root:

```bash
python3 guides/statistics/restore_project.py
```

This validates all SHA256SUMS files and snapshot structure without network access.

## Recreate a deleted Project

Set `GITHUB_PROJECT_TOKEN` in your local shell to a credential with Projects write access. Run the same command with `--apply`. The script creates a new Project and restores task associations, fields, values, statuses and visible-column views; existing issues and source commits are preserved. Local progress is saved for retries. Use a separate `--state` file to create another copy. A repository subset restores only the tasks in that subset; use the complete backup to restore the whole Project.

If issues were deleted too, set `GITHUB_REPOSITORY_TOKEN` to an issue-write credential and add `--recreate-missing-issues`. The tracker repository must exist. Recreated issues may receive new issue numbers; stable task IDs remain. Historical execution dates remain in fields, while GitHub creation/closure timestamps reflect restoration.

Repository associations are included in `project.json`. Set `GITHUB_REPOSITORY_TOKEN` (or a suitable `GITHUB_TOKEN`) during restoration to reconnect them automatically. Separate Project-write and repository credentials are supported. Comments, automations and permission grants are outside this reconstruction. Read [METHOD.md](METHOD.md) for estimates, unknown values and measurement boundaries. Git versions these files; no dated folder is needed.
