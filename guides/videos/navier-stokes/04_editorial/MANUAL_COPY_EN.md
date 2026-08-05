# Concise manual copy

Use this block as the introduction to a longer description or manual page.

## Why was this repository created?

The experiment asks whether an AI can navigate a mathematical research problem
outside its human operator's expertise while leaving enough formal checkpoints
for experts to inspect the work. The human supplied the objective, constraints
and continued direction. Codex selected and implemented the computational
route, wrote and ran the C++, interpreted failures, proposed narrower lemma
candidates, and documented the resulting proof frontier.

The main practical goal is time to first obstruction. A cheap scaling check,
identity test or small adversarial search may reject a weak lemma before a
researcher spends weeks or months trying to prove it. Empirical data may seed
the loop, but it is optional. The process starts with the cheapest reproducible
test and gives more compute only to claims that survive.

## What is the question?

The three-dimensional Navier–Stokes equations describe fluid motion. Vortices
can stretch structures into ever finer scales while viscosity smooths them.
The million-dollar question, in plain language, is whether a flow that starts
smooth must remain smooth forever—or whether an infinitely rough point can
form.

## Why does that matter outside mathematics?

Navier–Stokes equations are already used to model fluid behaviour in familiar
domains: atmospheric and weather-related flows, air around aircraft, and blood
flow in larger vessels. A solution to the Millennium problem would not
instantly produce better forecasts, aircraft or medical treatments. It would
settle a deeper foundation question: whether smooth three-dimensional
incompressible solutions are guaranteed to remain smooth, or whether the
mathematical model itself can develop a singularity. A proof gives certainty
and structure that experiments and finite simulations alone cannot provide.

## What did the two-day laboratory do?

It did not “test the whole problem.” It ran an AI-executed refinement loop:

1. turn the current proof gap into a candidate obligation;
2. collect code, certificates, traces and optional empirical artifacts;
3. attack the candidate with deterministic C++ and adversarial search;
4. let AI analyze the failure or surviving structure;
5. formulate a sharper candidate lemma and repeat.

The C++ computations provide reproducible gates. The AI connects one gate to
the next by reading the evidence and changing the mathematical target, code and
tests. A memory failure did not reset the campaign because saved states and
certificates preserved the frontier.

## What was the result?

Over 33 hours and 28 minutes, 16 lemma candidates were rejected and one partial
far-tail estimate was proved. The decisive cutoff-independent L4/PNT-12 step
remains open; surviving a finite computational stress test is not a proof. The
next PNT-12 decision cycle is estimated at roughly one to three days of focused
work. There is no honest estimate for a complete Clay solution.

That one-to-three-day figure is laptop pace, not a proof deadline. The campaign
used a 12-logical-CPU Intel Core i5-1335U. Ten comparable cluster nodes could
make independent cutoff/restart sweeps up to roughly ten times faster in their
ideal parallel portion. This is an engineering estimate rather than a measured
end-to-end benchmark; AI interpretation, memory-bound stages and mathematical
proof do not scale linearly.

This is a human-directed, AI-executed computational campaign, not an automated
proof. Experts are still needed to audit the mathematics and decide whether
the surviving structure matters.

## What work does it resemble?

Methodological kinship: Protas et al.'s adjoint-based search for extreme flows
and Tao's quantitative a priori bound viewpoint. In this laboratory,
adversarial search is used primarily to stress-test intermediate lemmas. This
is a comparison, not a claim that the work combines proofs by those authors.

## Why is there a one-million-dollar prize?

Clay originally selected seven Millennium Prize Problems and allocated US$1
million to each. The Poincaré Conjecture has been resolved, so six remain
unresolved. The award is not automatic if a claim sits unrefuted online. Before
CMI considers it, a proposed solution must be published in a Qualifying Outlet,
at least two years must pass, and it must gain general acceptance in the global
mathematics community. CMI does not accept direct solution submissions. After
the public qualification path, it conducts its own consideration.

## Required visual disclosure

> The fluid imagery is a concept visualization. The C++ laboratory actually
> produced code, saved states, certificates, and finite Fourier/Galerkin data.
