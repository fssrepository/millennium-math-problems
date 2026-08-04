# Millennium Math Problems

> **Human-directed, AI-executed computational proof search and
> lemma-falsification laboratory.** The current
> repository architecture, deterministic C++ laboratory, proof-state
> documentation, and public guides were produced under human direction by
> **OpenAI Codex Sol 5.6 at Extra High reasoning effort**. Replayable
> certificates, counterexamples, and restart states allow another AI or
> researcher to continue from the recorded proof frontier.

> **Ask the repository instead of reading it front to back.** Clone it, open
> the repository root in VS Code with the Codex plugin, and let Codex inspect
> the workspace. For example, ask: *"What does this project do? Explain the
> current proof route and PNT-12, and distinguish proved statements, finite
> numerical evidence, falsified routes, and open claims. Link every answer to
> the relevant source, certificate, state, or reproduction command."*

Public computational research laboratories organized by Millennium Prize
Problem. The active Navier–Stokes project uses deterministic C++ tests to reject
weak lemma routes early, shorten the research loop, and preserve useful
counterexamples as reproducible public knowledge.

> **Research status:** this repository has not solved the Navier–Stokes
> Millennium Problem. Finite numerical evidence is not a cutoff-independent
> proof.

## Public guides

| Guide | PDF | Slides | Purpose |
|---|---|---|---|
| Project overview | [Download PDF](https://raw.githubusercontent.com/fssrepository/millennium-math-problems/master/guides/pitches/output/Millennium_Math_Problems_Project_Overview_v1.0.0_EN.pdf) | [Download PPTX](https://raw.githubusercontent.com/fssrepository/millennium-math-problems/master/guides/pitches/output/Millennium_Math_Problems_Project_Overview_v1.0.0_EN.pptx) | What this repository does, the fast-falsification loop, optional empirical input, present evidence, and the proof roadmap. |
| Navier–Stokes brief | [Download PDF](https://raw.githubusercontent.com/fssrepository/millennium-math-problems/master/guides/problems/navier-stokes/output/Navier_Stokes_Brief_v1.0.0_EN.pdf) | [Download PPTX](https://raw.githubusercontent.com/fssrepository/millennium-math-problems/master/guides/problems/navier-stokes/output/Navier_Stokes_Brief_v1.0.0_EN.pptx) | The problem, why it is difficult, what the C++ laboratory tests, and how far the project has reached. |

## Problem portfolio

The definitions and solved/unsolved classification below link to the
[Clay Mathematics Institute](https://www.claymath.org/millennium-problems/).

| Problem | Goal | Status | Repository work |
|---|---|---|---|
| [Birch and Swinnerton–Dyer conjecture](https://www.claymath.org/millennium/birch-and-swinnerton-dyer-conjecture/) | Relate the rank of an elliptic curve to its L-function. | Unsolved | Planned. |
| [Hodge conjecture](https://www.claymath.org/millennium/hodge-conjecture/) | Characterize which topological classes arise from algebraic cycles. | Unsolved | Planned. |
| [Navier–Stokes existence and smoothness](https://www.claymath.org/millennium/navier-stokes-equation/) | Decide whether every smooth, physically reasonable 3D flow stays smooth, or exhibit finite-time breakdown. | **Unsolved · active here** | Built a spectral exact-gradient adversarial engine, proved a partial far-tail lemma, falsified standalone PNT-13 decorrelation, and is stress-testing PNT-12. <br> [`navier-stokes/`](navier-stokes/) |
| [P versus NP](https://www.claymath.org/millennium/p-vs-np/) | Decide whether every efficiently verifiable problem is efficiently solvable. | Unsolved | Planned. |
| [Poincaré conjecture](https://www.claymath.org/millennium/poincare-conjecture/) | Characterize the three-sphere among closed simply connected three-manifolds. | Solved | No active laboratory planned. |
| [Riemann hypothesis](https://www.claymath.org/millennium/riemann-hypothesis/) | Locate every nontrivial zero of the zeta function on the critical line. | Unsolved | Planned. |
| [Yang–Mills existence and mass gap](https://www.claymath.org/millennium/yang-mills-the-maths-gap/) | Construct the four-dimensional quantum theory and prove a positive mass gap. | Unsolved | Planned. |

## Research loop

| Stage | Activity | Persistent result |
|---|---|---|
| Direct | A scientist points Codex to a proof gap, the repository, and optional empirical or cluster data. | Concrete research target. |
| Formulate | Codex reads the proof state and proposes statements that would close the next gap. | Explicit candidate lemmas and assumptions. |
| Attack | Fast C++ gates check scaling and identities; exact-gradient adversaries search for hostile finite states. | Witnesses, traces, gradient checks, and certificates. |
| Falsify or refine | A counterexample is stored and used to narrow the next lemma; a survivor returns to analysis. | Failed-lemma record or sharper proof obligation. |
| Accumulate | New objectives, adversaries, importers, and restart points remain public. | Future work does not restart from zero. |

In the active campaign, this was one connected AI workflow rather than a fixed
simulation followed by manual analysis: Codex formulated candidates, wrote and
ran the C++, read failures, refined the next obligation, and documented the
frontier. The human operator set the objective and constraints, supplied the
hardware, reviewed checkpoints, and directed continuation. Experts are still
required to audit every mathematical claim.

The principal engineering metric is **time to first obstruction**: find the
cheapest reproducible reason a lemma cannot work before weeks or months are
spent trying to prove it. Optional empirical data can seed the same loop, but
the active Navier–Stokes campaign starts from the lowest open proof gap.
