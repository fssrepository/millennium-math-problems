# Two-entry local quartet

Combining the commutator identity with the projected normalization residual
reduces the complete local-local contribution to SLD-1P to two entries.

Let

```text
B=B_local(u,u),
S=<A u,B>,

K(u) = -<B,A B-B_local(u,A u)>
       +S^2/(2Z)
       +3S<A B,A u>/(2P),

G(u) = <A u,B_local(-B,u)>.
```

Here `K` contains the outer, advected, enstrophy-normalization, and
palinstrophy-normalization terms. `G` is the only remaining advecting Frechet
slot. The exact local part of the denominator-free numerator is

```text
N_local = 4 S^3 Z P (K(u)+G(u)).                     (LQR-4)
```

Thus a sufficient local sublemma for SLD-1P is the explicit two-entry bound

```text
4 S^3 Z P (K+G)
 <= A_local k0 (S^4 Z^2 P+B0 Z^3 P^4),              (SLD-1P-L)
```

with `A_local` independent of the Galerkin cutoff. The nonlocal and viscous
channels remain in their existing ledgers and cannot be silently charged to
this local statement.

On the optimized K6 state:

```text
K raw                              -0.195783730872
K normalized                       +7.4846687372e-4
G normalized                       -6.1297970097e-5
local normalized total             +6.8716890362e-4
N_local                            +0.00264627230688
independent polynomial N_local     +0.00264627230688
polynomial relative residual        3.20e-19
```

This is the current minimal local proof obligation. Applying a triangle
inequality to the original five entries is strictly weaker than estimating
`K+G` and is no longer the planned route.

[CLOSURE_TARGET.md](CLOSURE_TARGET.md) derives the exact sufficient estimate
for `K+G` and verifies its weighted Young reduction.
