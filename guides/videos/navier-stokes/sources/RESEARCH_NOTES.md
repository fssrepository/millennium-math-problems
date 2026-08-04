# Research and source notes

## What was inspected locally

### Session record

- Main campaign session: `019fbc49-c7f3-7051-8627-e920169afd1b`.
- Local archive:
  `/home/raxim/.codex/archived_sessions/rollout-2026-08-01T09-46-19-019fbc49-c7f3-7051-8627-e920169afd1b.jsonl`.
- Session message span: `2026-08-01T07:51:01Z` to
  `2026-08-03T00:10:11Z`. The result used by the film was reported at
  `2026-08-02T17:19:22Z`, exactly 33 h 28 min 21 s after the campaign began.
  See `../TIMELINE_EVIDENCE.md` for the event-level audit and screen mapping.
- The machine crash and resume were visible in the same record; the repository
  certificates and saved states allowed continuation instead of a restart from
  zero.

The session itself is private local provenance and is not copied into the
public bundle.

### Repository facts used in the film

- The problem and honest scope:
  `../../../navier-stokes/README.md`.
- Exact roadmap and still-open L4 step:
  `../../../navier-stokes/PROOF_PLAN.md`.
- Rejected candidates F000–F015:
  `../../../navier-stokes/proof/failed_lemmas.tsv`.
- Partial far-tail result and proof-state index:
  `../../../navier-stokes/proof/README.md` and the L4 lemma notes.
- Existing public three-page explainer:
  `../../problems/navier-stokes/output/Navier_Stokes_Brief_v1.0.0_EN.pdf`.
- Guide source text and original fluid hero:
  `../../problems/navier-stokes/build.py` and
  `../../problems/navier-stokes/assets/navier-stokes-flow.png`.

Snapshot checked on 2026-08-04:

- 301 tracked C/C++ source/header files under `navier-stokes/src`;
- 56,954 physical lines across those files by `wc -l`;
- 643 proof JSON files and 2,279 proof TSV files;
- 16 rejected candidates in `failed_lemmas.tsv`.

These scale numbers are useful for a description or making-of post, but were
left out of the 30-second narration to keep the scientific story legible.

## Public scientific source

- Clay Mathematics Institute, official Navier–Stokes problem description:
  https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf

## Protas / Tao comparison — verified scope

The comparison is useful, but it must be stated as **methodological kinship**.
It is not evidence that the local work literally combines their proofs.

### Protas side

- Kang, Yun and Protas optimize periodic 3D Navier–Stokes initial data with
  prescribed enstrophy to maximize finite-time enstrophy using large-scale
  adjoint-gradient methods; the computations find extreme trajectories but no
  singularity evidence: https://arxiv.org/abs/1909.00041
- Kang and Protas formulate an adjoint-based optimization around the critical
  Ladyzhenskaya–Prodi–Serrin quantity involving the time integral of the
  eighth power of the L4 norm: https://arxiv.org/abs/2110.06130
- A newer Ramírez–Protas study broadens the variational search across critical
  LPS criteria: https://arxiv.org/abs/2604.13338

This strongly resembles the local campaign's adversarial search for extreme
finite Fourier/Galerkin configurations and its use of exact gradients. Protas
optimizes flows/initial data; this campaign uses related adversarial machinery
primarily to attack intermediate lemma candidates.

### Tao side

- Tao's quantitative formulation relates qualitative periodic global
  regularity to the existence of quantitative local H1 a priori bounds:
  https://arxiv.org/abs/0710.1604
- Tao later derives explicit quantitative regularity information from control
  of a critical L3 norm: https://arxiv.org/abs/1908.04958

This resembles the architecture of seeking a numerical, cutoff-independent
critical-norm bound that would close a regularity route. It is a kinship of
proof strategy and quantitative viewpoint; the local repository does not
implement Tao's argument.

### Repository attribution audit

A name/citation search of the checked local source and guide tree found
Prodi–Serrin-related code and L4 objectives, but no explicit Protas or Tao
citation and no documented derivation from either body of work. Therefore the
safe public wording is:

> Módszertani rokonság: Protasék adjunkt-alapú szélsőáramlás-keresése és Tao
> kvantitatív a priori korlát-szemlélete. Itt az adverszárius keresés a köztes
> lemmákat stresszteszteli.

Avoid `Protas + Tao proof`, `built on their proof`, or any wording that implies
formal mathematical dependence not present in the repository.

## Platform sources checked on 2026-08-04

### Google Flow

- Model and feature matrix, including Veo 3.1 variants, 9:16 support, frames,
  ingredients and duration combinations:
  https://support.google.com/flow/answer/16352836?hl=en
- Official creation workflow for text, ingredients, start frames and model /
  aspect-ratio / duration selection:
  https://support.google.com/flow/answer/16353334?hl=en

### OpenArt

- Current video surface lists 9:16, Start + End Frame, Smart Shot, Kling 3.0,
  Veo 3.1 and Seedance, and describes camera prompts:
  https://openart.ai/ai-video-generator/
- Official Kling 3.0 Motion Control guide describing identity preservation and
  independent camera direction:
  https://openart.ai/blog/kling-3-motion-control-guide/
- Help center warning that model-generated text may be unreliable and noting
  standard clip lengths / external editing needs:
  https://openart.ai/help

### Higgsfield

- Official Kling 3.0 guide: 3–15 second clips, 2–6 scene structure,
  start/end-frame control, camera movement and macro strengths:
  https://higgsfield.ai/blog/Kling-3.0-is-on-Higgsfield-User-Guide-AI-Video-Generation
- Official 2026 comparison explaining Cinema Studio 3.5's genre, lighting,
  palette, camera, lens, focal-length and aperture controls, plus the platform's
  role guidance for Kling 3.0 / Veo 3.1 / Seedance 2.0:
  https://higgsfield.ai/blog/most-reliable-ai-video-generators-2026

## Editorial inference

The exact model recommendation is an editorial inference from those official
feature descriptions and the needs of this film. There is no universal “best”
video model. Here, the absence of people favors Veo for main plates; the fast
macro camera needs make Kling through Higgsfield the strongest retake path.
