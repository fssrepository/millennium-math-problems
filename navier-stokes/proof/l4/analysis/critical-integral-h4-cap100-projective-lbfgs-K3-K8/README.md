# Projected L-BFGS critical-integral continuation

This run starts from the optimized `K=3` projection of the earlier `K=8`
branch and successively lifts the state through `K=3,...,8`. At every cutoff it
performs 12 projected L-BFGS iterations on the exact discrete L4-A objective.
The initial energy is one, `nu=0.1`, `T=0.01`, `dt=0.002`, and the initial
homogeneous `H4^2` cap is 100.

| K | refined integral | increment | H4^2 | top-shell energy |
|---:|---:|---:|---:|---:|
| 3 | 6.118398892e-6 | — | 13.7629 | 8.1167e-5 |
| 4 | 6.119726945e-6 | 1.3281e-9 | 13.8715 | 3.4264e-7 |
| 5 | 6.119746587e-6 | 1.9642e-11 | 13.8758 | 1.4377e-9 |
| 6 | 6.119750514e-6 | 3.9265e-12 | 13.8809 | 2.4857e-12 |
| 7 | 6.119751027e-6 | 5.1311e-13 | 13.8808 | 5.0328e-14 |
| 8 | 6.119752080e-6 | 1.0528e-12 | 13.8794 | 1.4364e-13 |

The state-family certificate reports a maximum adjacent projection residual of
`2.2129e-4`. The branch is numerically cutoff-saturated and has a rapidly
decaying Fourier tail. This remains finite computational evidence: a proof
must derive a cutoff-independent estimate for every trajectory from the fixed
smooth datum, rather than assume persistence of this observed envelope.
