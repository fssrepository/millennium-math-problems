# Far-nonlocal adversary

This experiment splits the old nonlocal class into dyadic frequency gaps:

```text
gap 0: high/low <= 2                    (local)
gap 1: 2 < high/low <= 4                (near nonlocal)
gap 2+: high/low > 4                    (far nonlocal)
```

`critical-far-nonlocal-integral` maximizes

```text
J_far = integral_0^T |V_far|^4 / (Z P^3) dt
```

with the exact discrete RK4 adjoint. The nonzero K3 static gradient agrees with
central differences to `4.84e-13`. All tabulated values use `E(0)=1`,
`H4^2<=100`, `nu=0.1`, `T=0.01`, and a validated `dt=0.0005`/`0.00025` pair.

The result is in [`cutoff-summary.tsv`](cutoff-summary.tsv). The best K6 state
was projected separately to every lower cutoff so that optimizer progress is
not misread as cutoff growth:

| K | refined `J_far` | top-shell energy |
|---:|---:|---:|
| 3 | `2.11194e-11` | `1.05881e-3` |
| 4 | `3.20285e-11` | `2.68442e-5` |
| 5 | `3.21555e-11` | `6.81316e-8` |
| 6 | `3.21556e-11` | `1.37352e-10` |

The K5-to-K6 relative increment is about `1.33e-6`, while top-shell energy
drops by nearly three orders of magnitude. This optimized branch is consistent
with a smooth cutoff limit. It does not establish a uniform bound or rule out
another branch.

The dyadic ledger also identifies the hard interaction. On the separately
optimized full-nonlocal K5 state, gap 1 contributes signed stretching
`-0.371413`, whereas gaps 2+ contribute only `-8.92578e-5`. Within gap 1, the
low advecting and low advected roles dominate and the low target role partially
cancels them. This routes the proof decomposition to:

```text
L4.1a: bound the geometric tail gap >= 2;
L4.1b: combine gap 1 with the quasilocal estimate;
L4.2: control the resulting high/high/low transition and local block.
```

The split is exact, but these three bounds remain open.
