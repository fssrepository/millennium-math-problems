# Partitioned critical-integral adversaries

This experiment differentiates the local and nonlocal pieces of the proposed
L4-A density separately. For the fixed triad partition

```text
V = V_local + V_nonlocal,
local iff max(|k|,|p|,|q|) <= 2 min(|k|,|p|,|q|),
J_part = integral_0^T |V_part|^4 / (Z P^3) dt,
```

the discrete RK4 adjoint computes the exact gradient of `J_local` or
`J_nonlocal`. The deterministic central-difference test currently has relative
errors `4.23e-12` and `3.25e-12`, respectively. `TriadPartitioner` supplies the
same classification to the dynamics, trajectory diagnostics, and triad proof
checks.

## K=3 separation test

All runs use `E(0)=1`, `nu=0.1`, `T=0.01`, requested `dt=0.002`, and the
initial constraint `H4^2 <= 100`.

| optimized objective | refined `J_local` | refined `J_nonlocal` | final `H4^2` |
|---|---:|---:|---:|
| local | `4.38853e-6` | `7.20e-18` | `7.5668` |
| nonlocal, continued K3 | `1.22037e-10` | `1.61848e-6` | `100` |

The local optimizer selects a smooth, almost purely local state: `99.7355%`
of its energy is in shell one and only `4.98e-7` is in shell three. The
nonlocal optimizer moves `10.3055%` of the energy into shells two and three and
approaches the initial `H4` boundary. Thus the small nonlocal value observed on
the total-objective projective branch was branch-specific; it is not a
universal depletion law. The subsequent K4/K5 continuation is listed in
`../nonlocal-h4-cap100-converged/cutoff-summary.tsv`.

## Initial-H4 cap sweep

The preliminary replayable measurements are in
[`h4-cap-sweep.tsv`](h4-cap-sweep.tsv). Over caps `25, 50, 75, 100, 150, 200,
400`, the under-converged data gave the least-squares fit

```text
J_nonlocal ~ (H4^2 cap)^0.889697878
```

This exponent is rejected as an optimizer artifact. Extending the cap-100
search raised its objective from `2.10177e-7` to `1.61848e-6`, a factor of about
`7.7`, while the other cap points were not equivalently converged. No cap-law
exponent can be inferred from this sweep. The still-falsifiable analytical
candidate is

```text
L4.1-H candidate:
sup_N integral_0^T |V_N,nonlocal|^4 / (Z_N P_N^3) dt
    <= C(E(0), nu, T, ||u_0||_H4).
```

The right side may depend on the fixed smooth datum but not on Galerkin cutoff
`N`. A proof may use the initial frequency envelope and viscous smoothing; it
may not assume a uniform propagated `H4` solution bound, because that would be
circular. The next decisive machine check is a cutoff continuation of the
nonlocal objective at each fixed cap. The direct partitioned oracle is exact
but currently limits this search to small cutoffs; an FFT VJP for masked triad
classes is therefore the next performance kernel.

## Converged cutoff check at `H4^2 <= 100`

After removing optimizer progress by repeatedly projecting the best high-cutoff
state back to `K=3`, the continued branch gives

| cutoff | refined `J_nonlocal` | increment | top-shell energy |
|---:|---:|---:|---:|
| 3 | `1.61848084e-6` | - | `6.2453e-4` |
| 4 | `1.62065350e-6` | `2.1727e-9` | `1.0462e-6` |
| 5 | `1.62074608e-6` | `9.2579e-11` | `4.5120e-9` |

The K4-to-K5 relative increment is `5.71e-5`, and the top-shell energy drops
by more than two orders of magnitude per added shell. This branch is therefore
consistent with a smooth cutoff limit. It is one local optimum at one cap,
viscosity, and short horizon; it neither proves L4.1-H nor excludes a different
concentrating branch.

Source certificates and states:

- `proof/l4/adversary/partitioned-integral/l4-local-integral-h4-cap100-K3-lbfgs16.json`;
- `proof/l4/adversary/partitioned-integral/l4-nonlocal-integral-h4-cap*-K3-*.json`;
- `proof/l4/states/local-integral-h4-cap100-K3-lbfgs16/`;
- `proof/l4/states/nonlocal-integral-h4-cap-sweep/`;
- `local-K3.json` and `nonlocal-K3.json` in this directory.
