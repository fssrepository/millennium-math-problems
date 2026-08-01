# Transition-block obstruction

The moving gap

```text
m(Z) = m0 + ceil(log2(max(1,Z)))
```

reduces the geometric far-tail remainder to a linear `Z` term. It leaves at
most `O(log(1+Z))` local and transition dyadic bands. Merely counting these
bands does not close the estimate.

The undepleted local estimate has the form

```text
|V_local+transition|
    <= C log(1+Z) Z^(3/4) P^(3/4).
```

Young's inequality with conjugates `4/3` and `4` gives

```text
|V_local+transition|
    <= (nu/4) P + C nu^(-3) log(1+Z)^4 Z^3.
```

The energy identity controls `integral Z dt`. The normalized remainder relative
to that controlled term is

```text
Z^2 log(1+Z)^4,
```

which diverges as `Z` grows. A logarithmic number of bands changes no
polynomial exponent. `TransitionBlockScaling` verifies the exact rational
exponents and records the required pointwise depletion power `Z^(-1/2)`.

Therefore the transition/local block needs a genuinely dynamical or geometric
depletion estimate. Band counting alone is rejected as candidate F005.
