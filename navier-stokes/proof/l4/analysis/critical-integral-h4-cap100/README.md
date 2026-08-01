# Critical-integral H4 continuation

This experiment directly maximizes the discrete L4-A candidate

```text
J_K(u_0) = integral_0^0.01 D_K(t)^4 Z_K(t)^2 dt
```

at unit initial energy, viscosity `0.1`, and initial homogeneous `H4^2 <= 100`.
The optimizer uses the exact discrete RK4/trapezoid adjoint, fixed-energy
Riemannian projection, Sobolev tangent projection, and backtracking retraction.

## Result

| K | search start | search finish | refined integral | initial H4^2 |
|---:|---:|---:|---:|---:|
| 1 | 6.6328e-9 | 3.6579e-6 | 3.6579e-6 | 3.6856 |
| 2 | 3.6938e-6 | 5.8919e-6 | 5.8919e-6 | 11.0343 |
| 3 | 5.8958e-6 | 5.9996e-6 | 5.9996e-6 | 12.9023 |
| 4 | 5.9997e-6 | 6.0060e-6 | 6.0060e-6 | 13.2261 |
| 5 | 6.0060e-6 | 6.0083e-6 | 6.0084e-6 | 13.3523 |
| 6 | 6.0083e-6 | 6.0095e-6 | 6.0095e-6 | 13.3868 |
| 7 | 6.0095e-6 | 6.0100e-6 | 6.0101e-6 | 13.3856 |
| 8 | 6.0100e-6 | 6.0104e-6 | 6.0104e-6 | 13.2794 |

A separate 24-step warm continuation at `K=8` increases the coarse search
objective from `6.01043e-6` to `6.04548e-6`; half-step replay gives
`6.04549e-6`, a relative time-discretization difference of `1.18e-6`. Its
energy-balance residual is `-8.29e-12` and `H4^2=15.7589`.

At relative energy threshold `1e-8`, the refined state has 294 active modes and
highest active shell three. Shell one carries `98.93%` of the energy. The
cutoff shell carries `1.49e-11`, so this branch is not using the `K=8` boundary
to inflate the objective.

## Logical status

The observed branch is consistent with a smooth cutoff limit on this short
trajectory. It is neither a proof of a cutoff-uniform bound nor evidence for a
singularity. The remaining test parameters are longer time, lower viscosity,
and independent multistarts. A proof still needs an analytical frequency
envelope controlled from the fixed smooth initial datum without assuming the
high-norm bound that L4 is intended to establish.
