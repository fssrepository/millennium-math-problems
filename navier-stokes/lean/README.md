# MATH-19: Lean PNT audit

Lean proves that the unrestricted PNT-2/7/8/12/14 upper bounds are false on an
explicit six-mode Fourier family. The full cutoff-one output cube, Cauchy
factorization, finite Gram-Schur estimate and concrete Gram row are checked.
These are finite mathematical counterexamples, not a result about global
Navier–Stokes regularity.

## Run

From the repository root:

```bash
./navier-stokes/lean/verify.sh
```

Requirements: Linux, elan, Git, Python 3, Bash, and a C++20 compiler with OpenMP
and GNU-compatible section garbage collection (tested with GCC). The first run
downloads the pinned Lean toolchain, dependencies and the required mathlib cache.
Later runs reuse them. The existing C++ application build remains independent
of Lean. `LAKE` and `CXX` may select the Lake and C++ executables.

If elan is not installed:

```bash
curl --fail --silent --show-error --location --retry 3 \
  https://raw.githubusercontent.com/leanprover/elan/e9268f8b213e09a99768befe548049c046074225/elan-init.sh \
  -o /tmp/math19-elan-init.sh
sh /tmp/math19-elan-init.sh -y --no-modify-path --default-toolchain none
```

The verifier also finds elan at `$HOME/.elan/bin/lake`; no shell-profile or
global default-toolchain change is required. Lean is pinned to 4.35.0. Mathlib
is pinned to commit `56a005855ea17e0b0958f73ec801ae140b0a1679` (v4.35.0);
`lake-manifest.json` pins the transitive dependencies. Verification does not
run `lake update`.

## Folder structure

```text
lean/
├── NavierStokesAudit.lean           public proof-library imports
├── NavierStokesAudit/
│   ├── Algebra/Normalization.lean   scalar formulas and scaling laws
│   ├── Fourier/
│   │   ├── SixMode.lean             original admissible PNT-2 witness
│   │   ├── Cube.lean                full output cube and height shells
│   │   └── Embedding.lean           zero extension, full norms and pairings
│   └── PNT/
│       ├── Cauchy.lean              alignment and squared majorant
│       ├── Schur.lean               general finite Gram-Schur inequality
│       └── Majorants.lean           concrete PNT-7/8/12/14 counterexamples
├── ProofAudit.lean                 public audit-command import
├── ProofAudit/Commands.lean        transitive axiom checks
├── tests/
│   ├── lean/Audit.lean             required results and axiom audit
│   ├── lean/rejections/            deliberately rejected proofs
│   └── cpp/witness.cpp             comparisons with original C++ owners
├── scripts/check-cpp.sh            focused C++ build and comparison
├── verify.sh                      single verification entry point
├── IMPLEMENTATION_PLAN.md          ordered, self-contained completion record
├── README.md
├── lean-toolchain
├── lakefile.toml
└── lake-manifest.json
```

Generated dependencies and build products stay in `.lake`; verification logs
and the C++ executable stay in `.verification`. Both directories are ignored.
Existing theorem namespaces remain stable after moving their source files.

## What is proved

| Result | Scope |
| --- | --- |
| Core-tail factorization | The PNT-1 three-term and two-term identities over the real numbers. |
| Open quotient scaling | PNT-2 scales as `1/a`, with the explicit amplitude action and nonzero denominators. |
| Squared-majorant scaling | PNT-7/Schur scales as `1/a^2` when the Gram row is amplitude-invariant. |
| Fixed-energy restriction | For positive energy and amplitude, preserving energy forces `a=1`. |
| Concrete PNT-2 counterexample | Admissibility, geometry, convolution, denominators and exact quotient are proved. |
| Cauchy chain | Finite Cauchy-Schwarz; alignment ≤ 1 given the two pairing bounds; selected quotient = majorant × sqrt(alignment); the exact squared-majorant formula. Positivity assumptions appear in each theorem. |
| Finite Gram-Schur estimate (PNT-11) | For any finite real vector family, the squared norm of the sum is bounded by a normalized absolute Gram-row upper bound times the diagonal sum. Zero vectors and zero rows are handled. |
| Concrete PNT-7/8/12/14 counterexamples | The actual full-output norms, shell decomposition, diagonal sum and maximum/fixed Gram row yield a positive majorant. At H=1 every finite bound is exceeded on the admissible amplitude ray. |

The majorant counterexamples no longer assume an unrealized positive base value:
Lean computes it from the Fourier coefficients. The H=1 slice is sufficient to
refute a bound claimed uniformly over all allowed H≥1 and cutoffs. No assertion
is made about a fixed-energy supremum. This package does not prove C++ program
or gradient correctness, the separate PNT-13 height-decorrelation claim, the
positive local-density-derivative candidate, or continuous-PDE regularity.

## The exact witness

The real parts of all Fourier coefficients are zero. At positive amplitude
`a`, the imaginary parts on the positive representatives are:

| Wave | Imaginary coefficient |
| --- | --- |
| `(1,0,0)` | `a (0,1,0)` |
| `(0,1,1)` | `a (1,0,0)` |
| `(1,1,1)` | `a (1,-1,0)` |

At each opposite wave the imaginary coefficient is negated. All other
coefficients, including the zero mode, vanish. The Lean model proves that the
six waves are distinct and nonzero, their opposite-wave map is correct, and
the coefficients obey reality and divergence-free constraints.

The six-mode module computes the Fourier convolution and Leray projection on
the velocity support for PNT-2. The cube module additionally computes all 27
cutoff-one outputs, with the zero output masked. The zero extension is proved
mean-zero, real and divergence-free; its norms and pairings agree with the
compressed input. In addition to the six supported outputs, the waves
`±(1,-1,-1)` contribute to the advection norms. Omitting these would give the
wrong majorant. Every retained input interaction is local and has primitive
squared-length height 3. The height-one core is empty and tail shell 2 is the
only nonzero shell; all other natural-number shell indices are proved zero.
The normalization agrees with the laboratory's full
Fourier sum convention, with no extra factor 1/2 in energy.

```text
E = 8 a^2       Z = 18 a^2       P = 46 a^2
S_full = S_tail = 2 a^3          T_tail = 6 a^3
S_core = T_core = 0
PNT-2 = 3 |S_full| |S_tail T_tail| / (2 Z^2 P^3)
      = 1 / (876024 a)
```

`NavierStokesAudit.FourierWitness.exceeds_every_bound` proves that every real
candidate bound is exceeded by an admissible positive-amplitude member of this
family. `NavierStokesAudit.FourierWitness.no_uniform_bound` consequently refutes
a uniform bound over all admissible states on this support with positive
enstrophy and palinstrophy. No initial energy is fixed in that quantifier.

For the full cutoff-one output cube, Lean additionally proves:

```text
B1 = 22 a^4       T2 = D_H = 50 a^4
Gram rows (shells 0,1,2) = (0,0,1) for a > 0
R_H = R_(H,2) = 1        sqrt(H) = 1 at H = 1
PNT-7 = PNT-8(H=1) = PNT-12(H=1) = PNT-14(H=1,row=2)
      = 275 / (33366002112 a^2)
```

`NavierStokesAudit.Majorants.exceeds_every_bound` supplies one admissible
amplitude exceeding any proposed real bound in all four formulas simultaneously.
The `no_cauchy_bound`, `no_schur_bound` and `no_fixed_row_bound` theorems state the
corresponding negated uniform bounds. The Gram row is evaluated from its actual
normalized inner products, including zero rows, rather than assigned as a premise.

## Verification and trust boundary

The single command performs these checks in order:

1. Every proof module is imported by the audited root.
2. Lake builds the proof library and the audit command.
3. Lean's transitive axiom collector checks every declaration in the project
   namespace, including private declarations. Only `propext`,
   `Classical.choice` and `Quot.sound` are allowed. The collector uses Lean's
   exported axiom metadata, which accounts for erased imported proof bodies.
4. Lean rejects false open and squared-majorant bounds at amplitude 1/2. It also
   rejects a deliberate unfinished proof and an injected axiom. These four
   fixtures under `tests/lean/rejections`
   intentionally fail and are excluded from the proof library.
5. A small C++ driver calls the existing spectral and normalization owners for
   amplitudes 1, 1/2 and 1/4 in complete cutoff-one and cutoff-two cubes. It
   checks the pairings, quotient and existing amplitude/scope regression. In the
   cutoff-one cube it also calls the original Cauchy and Schur owners and checks
   B1, T2, D_H, R_H and both majorant formulas at all three amplitudes.

The real-number algebra and finite computations produce proof terms checked
by Lean; no floating-point certificate or `native_decide` is used as a premise.
The usual Lean kernel, imported library artifacts and toolchain remain the
trusted software base. This is not an independently checked second-kernel
certificate or a formal equivalence proof between Lean and compiled C++.

The C++ comparison is a numerical consistency regression using the original
implementation, not a mathematical premise. Its tolerance is relative
`1e-14`, with a `1e-20` scale floor for expected zero values. The full C++
research self-test is not rerun by this focused verifier. The C++ evaluators
have floating-point small-denominator safeguards; the finite comparison samples
stay above those thresholds. The theorem for arbitrarily small positive real
amplitudes is a Lean result and does not assert numerical fidelity below the
C++ safeguards.

Successful output ends with:

```text
MATH-19 PASS: unrestricted PNT-2/7/8/12/14 bounds disproved; proof audit and C++ consistency checks passed.
```

Logs and the C++ probe executable are generated under `.verification`, which
is ignored by Git. The GitHub Actions workflow `Lean PNT audit` invokes the
same command. Local verification passed; a hosted Actions run is not claimed
by the local completion record.

The ordered implementation and its completion results are recorded in
`IMPLEMENTATION_PLAN.md`, a self-contained numbered document.
