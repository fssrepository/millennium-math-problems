# MATH-19: Lean introduction, ordered implementation and completion record

Scope: add a reproducible Lean 4 check of the withdrawn PNT amplitude claims.
The first unconditional counterexample targets PNT-2. This is not a proof or
disproof of global Navier–Stokes regularity, nor an audit of every repository
lemma. No numerical search campaign is part of this task.

Work proceeds in the numbered order below. Finish and verify each prerequisite
before its consumers. Completed stages are not restarted. This document contains
the complete plan and results without cross-links or a separate task ledger.

1. [x] Assess the existing definitions and fix the scope.
   The C++ normalization objective uses
   `3 |S_full| |s_c t_t + s_t t_c + s_t t_t| / (2 Z^2 P^3)`.
   Its squared Cauchy/Schur majorants have numerator amplitude degree 14 and
   denominator degree 16. The open objective has degrees 9 and 10.
   Positive amplitude rescaling preserves finite support, Fourier reality,
   incompressibility, locality and primitive projective shape. Fixed-energy
   normalization excludes this amplitude ray. There is no existing Lean
   project, CI workflow or installed Lean toolchain to extend. The MyScoutee
   application template contains no equivalent mathematical implementation.
   Acceptance: actual formulas and admissibility conditions are identified;
   no application/template change is needed for this independent research repo.

2. [x] Establish the reproducible Lean environment.
   Pin Lean 4.35.0 and the matching mathlib release/commit; save the dependency
   manifest; ignore generated Lake files. Install elan without changing the
   user's shell profiles or global default toolchain. Verify a minimal import
   before implementing mathematical modules.
   Result: Lean 4.35.0 and mathlib commit
   `56a005855ea17e0b0958f73ec801ae140b0a1679` installed; the minimal real-number
   import builds without warnings. Dependency versions are saved in the
   manifest. Only the required mathlib cache modules are needed. The initial
   GitHub HTTP/2 fetch failure was resolved with a shallow HTTP/1.1 fetch.

3. [x] Formalize the scalar PNT definitions and amplitude argument.
   Prove the two core-tail factorizations, degree-minus-one scaling for the
   open objective and degree-minus-two scaling for the squared majorant.
   Prove that a positive value with either scaling law defeats every finite
   uniform bound on a positive amplitude ray. State every nonzero/positivity
   assumption explicitly. Prove the fixed-energy restriction separately.
   Acceptance: Lean checks these results, with no unfinished proofs or new
   axioms. These are general algebraic implications, not yet a Fourier witness.
   Result: the scalar module passes Lean with warnings treated as errors.
   Both scaling laws, both factorization identities, explicit amplitudes
   defeating any proposed bound, and the fixed-energy restriction are proved.

4. [x] Supply a concrete admissible finite Fourier witness for PNT-2.
   Use a small exact, divergence-free, conjugate-symmetric spectrum, verify
   its local/core-tail classification and evaluate its pairings in Lean.
   Connect the computed positive objective to the amplitude theorem. The final
   theorem must refute a uniform bound on this explicit family, without assuming
   the desired nonzero numerator or importing a floating-point certificate.
   Compare the same witness with the existing C++ objective as a separate
   implementation consistency check; C++ output is not a proof premise.
   Result: Lean verifies the six modes, conjugate symmetry, divergence freedom,
   locality, primitive height 3, empty height-one core, and the projected
   convolution. It proves `E=8a^2`, `Z=18a^2`, `P=46a^2`, `S=2a^3`,
   `T=6a^3`, and `PNT-2=1/(876024a)` for every positive `a`, then proves
   that no finite uniform bound exists. The existing C++ objective agrees for
   amplitudes 1, 1/2 and 1/4 in complete cutoff-one and cutoff-two cubes.

5. [x] Integrate a single verification command and continuous integration.
   Build all proof modules; inspect the axioms of the public results; reject
   unfinished proofs and extra trusted axioms. Include a deliberately false
   proof as an expected-failure regression so acceptance and rejection are both
   exercised. Reuse the same command locally and in GitHub Actions. Preserve
   the independent, dependency-free C++ build.
   Result: `verify.sh` passes the proof build, the 123-declaration transitive
   axiom audit, all three expected rejections (false bound, unfinished proof,
   injected axiom), six C++ witness comparisons, and the C++ amplitude/scope
   regression. The GitHub Actions workflow uses this exact command with pinned
   action revisions. A hosted workflow run is not part of the local evidence.

6. [x] Run final checks and finish the usage/results documentation.
   Record the exact commands, theorem statements, scope, trust boundary and
   observed outcomes here and in the Lean README. Report which PNT statement is
   disproved, which results are conditional, and what is outside this task.
   Final acceptance: the proof build and audit pass, the false proof is rejected,
   the existing C++ witness agrees, and every step in this plan is complete.
   Result: the complete local verifier passed on 2026-10-25. The two proof
   modules contain 27 named theorems; the audit covers 123 declarations,
   including generated dependencies. Shell syntax, workflow structure and all
   manifest commit pins were checked. The README records installation,
   reproduction, formulas, trust assumptions, expected failures and exact scope.
   The canonical PNT note records the new formal counterexample. The complete
   PNT-7/8/12/14 analytical chain and C++ correctness are explicitly outside the
   completed result. Task and commit use the identifier MATH-19.

Continuation requested 2026-10-25, started 16:25 UTC. Completed stages 1–6 remain
closed. The continuation proves the remaining unrestricted majorant claims false
on an actual Fourier witness and formalizes the finite Schur inequality; it does
not assert correctness of the optimizer, gradients, or the full PDE.

7. [x] Organize the existing library by responsibility, preserving its proofs.
   Algebra contains normalization identities; Fourier contains spectral models;
   PNT contains the majorant conclusions. Put proof auditing in ProofAudit,
   regression fixtures in tests/lean and tests/cpp, helpers in scripts.
   Keep verify.sh as the single public entry point.
   Result: existing proofs rebuilt successfully in the new structure. No
   mathematical statement or namespace was changed by the move.

8. [x] Extend the Fourier model to every output in the cutoff-one cube.
   Compute the full projected convolution, including outputs outside velocity
   support. Prove the selected H1 norm and tail H2 norm exactly, and classify
   every contributing interaction by primitive height and dyadic shell.
   Result: Lean proves cube coverage, support embedding, full admissibility and
   full-cube pairings, B1=22a^4 and T2=50a^4. Two advection outputs outside the
   velocity support are included. Every interaction has primitive height 3;
   shell 2 is the only possible nonzero shell, including among all natural
   shell indices. The first full-output proofs passed about 8 minutes after
   continuation started.

9. [x] Prove finite Gram-Schur control and the concrete PNT-7/8/12/14 conclusions.
   Handle zero rows explicitly. Evaluate the concrete Gram row and diagonal sum;
   prove positivity and an admissible witness exceeding every proposed bound.
   Result: the Cauchy alignment identity and squared majorant, the general
   finite Gram-Schur vector estimate (including zero rows), the actual row
   R=1 and D=T2=50a^4 are proved. At H=1 all four majorants equal
   275/(33366002112a^2). Their common positive Fourier witness exceeds every
   real bound. All proof modules built successfully by approximately 16:40 UTC.
   PNT-13 standalone height-shell decorrelation, gradients and continuous-PDE
   conclusions are separate claims and are not asserted by these theorems.

10. [x] Extend the single verifier and close the continuation.
   Audit all new proofs, reject a false majorant bound, compare exact constants
   with the existing C++ Cauchy/Schur owners, update documentation and MATH-19,
   and create the MATH-19 commit. Record measured times and precise scope.
   Result: the complete verifier passed on 2026-10-25 by 16:41 UTC, about
   16 minutes after the continuation began. It checks 75 named theorems in
   7 modules and audits 259 declarations with only propext, Classical.choice
   and Quot.sound allowed. All four negative fixtures were rejected for the
   intended reason. The original C++ owners agree at amplitudes 1, 1/2 and
   1/4, including B1, T2, D_H, R_H and the squared majorants at cutoff one.
   The earlier 2–3 day / multiweek estimates were not supported by observed
   work; this continuation's measured implementation-and-verification time
   is about 16 minutes, with documentation/task/commit closure afterwards.
   Folder ownership and trust boundaries are documented. The commit remains
   local; no hosted Actions run is claimed. MATH-19 is the task/commit identifier.
