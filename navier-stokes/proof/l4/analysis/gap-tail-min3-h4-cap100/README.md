# Dyadic gap-three tail adversary

This experiment generalizes the exact triad mask to a selectable tail:

```text
gap >= m  iff  largest |k| / smallest |k| > 2^m.
```

The command-line objective `critical-gap-tail-integral` maximizes

```text
J_m = integral_0^T |V_gap>=m|^4 / (Z P^3) dt
```

with the discrete RK4 adjoint. The run uses `m=3`, `E(0)=1`, initial
`H4^2<=100`, `nu=0.1`, `T=0.01`, 12 CPU workers, and projected L-BFGS. Every
accepted result below was repeated with `dt=0.0005` and `dt/2`; the relative
time-step errors are below `6.5e-5`.

The first eight steps exposed a scale bug in the old optimizer tolerances:
they improved the objective from `1.25e-43` to `1.32e-21`. Relative
machine-epsilon acceptance replaced the absolute small-objective cutoff. Three
four-step checkpoint blocks then reached `6.07829e-21` at K6. Four K7 steps
reached `1.29741e-20`.

Two controls separate cutoff effects from optimizer progress:

| source state | measured K | refined `J_gap>=3` | retention |
|---|---:|---:|---:|
| K6 winner | 5 projection | `1.49098e-22` | `2.45%` of K6 |
| K6 winner | 6 | `6.07829e-21` | `100%` |
| K6 winner | 7 zero lift | `6.00481e-21` | `98.79%` of K6 |
| K7 winner | 6 projection | `1.05251e-20` | `81.12%` of K7 |
| K7 winner | 7 | `1.29741e-20` | `100%` |
| K7 winner | 8 zero lift | `1.29712e-20` | `99.9775%` of K7 |

K5 barely contains any gap-three triads, which explains its poor projection.
Once that band opens, the K6-to-K8 controls are stable. The K7 state has
top-shell energy `1.67e-8`; its dyadic ledger records signed gap-three
stretching `5.6910e-5` and absolute pair stretching `3.8916e-4`.

The exact low-advecting commutator now isolates an analytical gain. With
`k=p+q` and

```text
T(p,q,k) = Re <u_k, i(q dot u_p) u_q>,
```

reality and `p dot u_p=0` give

```text
T(p,q,k) + T(p,-k,-q) = 0,
|k|^2 T(p,q,k) + |q|^2 T(p,-k,-q)
    = (|k|^2-|q|^2) T(p,q,k).
```

Since `||k|^2-|q|^2| <= |p|(|k|+|q|)`, a gap-`m` pair gains one factor
comparable to `2^-m` relative to the two-derivative unpaired weight. On the K7
gap-three state, pairing reduces the absolute low-advecting contribution from
`3.53962e-4` to `2.58109e-5`, a ratio `0.07292`. The weighted identity residual
is `2.05e-20`, and its maximum normalized frequency inequality ratio is
`0.8512`. This is exact finite Fourier algebra, not an empirical fit.

`TriadTailEnvelope` now covers all three placements of the unique low wave.
For the K7 gap-three state its exact finite-mode ledger is:

| low-wave role | signed stretching | absolute stretching | coefficient-free envelope | max amplitude ratio | max frequency ratio |
|---|---:|---:|---:|---:|---:|
| advecting | `2.36079e-5` | `2.58109e-5` | `7.82340e-5` | `0.98347` | `0.80623` |
| advected | `3.38163e-5` | `3.46702e-5` | `1.39892e-4` | `0.94175` | `1.00000` |
| target | `-5.14118e-7` | `5.26813e-7` | `1.72405e-5` | `0.11416` | `1.00000` |

The three signed rows sum to the complete gap-three stretching. The finite
certificate verifies the amplitude and frequency inequalities term by term;
it does not supply the cutoff-independent shell summation.

For this state `Z=1.14701`. With base gap two, the dynamic rule

```text
m(Z) = 2 + ceil(log2(max(1,Z)))
```

selects gap three and gives
`2^(-2m) Z^3 / Z = 0.0205569 <= 2^(-4)`. This is one state-level check of the
algebraic moving-gap inequality, not a trajectory regularity assumption.

For comparison, the separately optimized `gap>=2` integral is
`3.21556e-11`. The current gap-three value is smaller by a factor of about
`2.48e9`. This is evidence of rapid tail decay on the found branches, not a
proof of an exponent: signed cancellation is raised to the fourth power, the
search is not globally certified, and only short-time finite cutoffs were
tested.

The proof target is still a cutoff-independent summable paraproduct envelope
whose constants use only already controlled quantities and the fixed smooth
initial datum. The numerical branch decides which estimate to attempt; it
does not replace that estimate.
