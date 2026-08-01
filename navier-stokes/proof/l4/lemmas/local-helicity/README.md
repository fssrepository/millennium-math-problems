# Local helical decomposition

`HelicalTriadLedger` decomposes every divergence-free Fourier coefficient into
curl eigenvectors

```text
i k cross h_s(k) = s |k| h_s(k),   s in {-1,+1},
```

and reconstructs vortex stretching from all eight sign sectors
`(s_p,s_q,s_k)`. The decomposition is exact: the current certificate has
maximum velocity, total-stretching, and local-stretching reconstruction
residuals below `2.5e-20`.

## Adversarial-state observation

For the K3 local-objective endpoint,

```text
V_local              = -2.2778266270e-1
V_local,homochiral   = -3.7146823897e-6
V_local,heterochiral = -2.2777894802e-1
E_plus / E_minus     = 0.4999997884 / 0.5000002116
```

Thus the signed local extremizer is almost entirely heterochiral. Its
homochiral absolute interaction sum is `6.46194e-1`, so the small signed value
comes from cancellation rather than absent homochiral interactions. The
separately optimized nonlocal endpoint has the same qualitative local split:
`V_local,homochiral=-5.18518e-4` and
`V_local,heterochiral=-3.41206e-2`.

## Rejected identity

The observed cancellation is not universal. Projecting deterministic random
divergence-free states onto pure positive or pure negative helicity leaves no
heterochiral sectors but produces nonzero local enstrophy stretching. The
certificate reaches

```text
max |V_local,homochiral| = 6.35947e-3.
```

Therefore candidate F006, exact homochiral cancellation, is rejected. The
surviving analytical direction is narrower: identify a quantitative dynamic
depletion of the heterochiral local sectors, or reject that proposal with a
sector-targeted adversary. The finite endpoint pattern is routing evidence,
not a regularity proof.
