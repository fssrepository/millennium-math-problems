# Effective-signature adversary

This search tests whether the conditional effective-count reduction LSF-3 can
be closed by an algebraic estimate valid for every fixed-energy Fourier state.
It does not optimize the Navier--Stokes trajectory objective. Its sole purpose
is to reject an overstrong lemma before analytical work is built on it.

For each cutoff it generates three divergence-free, real-valued spectral
profiles:

- the existing algebraically decaying random state;
- a flat-spectrum state obtained by exactly undoing that generator's decay;
- a flat state supported in the outer half of the Fourier cube.

Every state is normalized to unit energy. `LocalTriadSymmetrizer` computes the
complete-triad signed transfers `V_sigma` and

```text
N_eff = (sum_sigma |V_sigma|)^2 / sum_sigma |V_sigma|^2.
```

The 24-sample, 12-worker K1--K6 run reports the following maxima:

```text
K   decaying       flat       outer-half-flat
1     1.98495      1.99883          1.99938
2    15.6244      16.6872           7.96149
3    60.9765      67.4531          60.2996
4   133.007      245.076          159.764
5   312.441      657.398          577.079
6   610.246     1748.62          1228.99
```

The flat and outer-half-flat finite-range log slopes are `4.223` and `4.554`.
Finite sampling does not establish an asymptotic lower bound, but these values
decisively fail the project's adversarial acceptance test for an unrestricted
sublinear-count lemma. In particular, the broad optimized K3 endpoint's
`N_eff=1.1902` is special trajectory structure, not an identity following
from incompressibility and fixed energy alone.

The next viable candidate must use time integration, viscous weighting,
smooth initial-data persistence, or an amplitude-weighted signature measure.
It cannot assert a pointwise sublinear `N_eff` bound for arbitrary Fourier
states.

The tighter signed quantity is

```text
A_sig = |sum_sigma V_sigma| / (sum_sigma |V_sigma|^2)^(1/2).
```

Unlike `sqrt(N_eff)`, it preserves cancellation between signatures. At K6 the
largest sampled `A_sig` values are only `2.33549` (flat) and `2.73471`
(outer-half-flat). Their K2--K6 fitted slopes are `0.48310` and `0.25607`.
Both lie below the closure threshold `1/2`; this survives the finite search
but remains an open analytical candidate.

The machine-readable result is
[`../local-signature-adversary-K1-K6.json`](../local-signature-adversary-K1-K6.json).
Reproduce it with

```bash
./build/navier_stokes_lab local-signature-adversary \
  --min-cutoff 1 --max-cutoff 6 --samples 24 --workers 12 \
  --certificate proof/l4/analysis/local-signature-adversary-K1-K6.json
```
