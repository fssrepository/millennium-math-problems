# Converged nonlocal cutoff branch

This directory separates optimizer progress from cutoff growth for the exact
nonlocal L4-A objective at `E(0)=1`, `H4^2<=100`, `nu=0.1`, and `T=0.01`.
High-cutoff winners were projected back to `K=3` until the low-mode branch was
nearly stationary, then lifted through `K=4` and `K=5`.

The machine-readable result is [`cutoff-summary.tsv`](cutoff-summary.tsv).
The refined objective increments fall from `2.1727e-9` at K4 to `9.2579e-11`
at K5, while top-shell energy falls from `1.0462e-6` to `4.5120e-9`. The
observed branch is consistent with a smooth cutoff limit. It is not a global
optimization certificate and does not prove a cutoff-independent bound.

`K3.json`, `K4.json`, and `K5.json` contain shell and Sobolev diagnostics. The
source optimizer certificates are:

- `proof/l4/adversary/l4-nonlocal-integral-h4-cap100-K3-convergence128.json`;
- `proof/l4/adversary/l4-nonlocal-integral-h4-cap100-K3-to-K4-test.json`;
- `proof/l4/adversary/l4-nonlocal-integral-h4-cap100-K4-to-K5-test.json`.
