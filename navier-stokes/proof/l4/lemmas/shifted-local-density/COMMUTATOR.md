# Outer--advected commutator

Let `B=B_local(u,u)` and `A=-Delta`. In the local nonlinear direction `h=-B`,
the outer and advected Frechet slots are

```text
O = -<A B,B>,
D = -<A u,B_local(u,B)>.
```

The local triad mask is symmetric under the transport pairing. Periodicity and
`div u=0` therefore preserve the usual skew identity

```text
<f,B_local(u,g)> = -<B_local(u,f),g>.
```

Taking `f=A u` and `g=B` gives the exact combined identity

```text
O+D
 = -<B, A B-B_local(u,A u)>.                         (LQC-1)
```

In Fourier variables, the multiplier in the commutator is `|k|^2-|q|^2`,
where `p+q=k` and `p` is the advecting frequency. It satisfies

```text
abs(|k|^2-|q|^2)
 = abs(p dot (k+q))
 <= |p|(|k|+|q|).                                    (LQC-2)
```

`LocalQuarticCommutator` evaluates both sides of LQC-1 independently and
checks LQC-2 on every included interaction. Equality in LQC-2 is possible, so
the executable allows only a `1e-15` floating-point slack around one.

On the optimized K6 state:

```text
outer slot                         -0.603413126167
advected slot                      +0.325899010744
combined                           -0.277514115423
negative commutator pairing        -0.277514115423
relative identity error             1.95e-19
Cauchy ratio                         0.860459094858
maximum normalized symbol ratio      1
```

LQC-1 removes one independent large term from the proof ledger. It does not
by itself determine the sign after multiplication by `4S^3`, whose sign
depends on the signed local stretching `S`.
