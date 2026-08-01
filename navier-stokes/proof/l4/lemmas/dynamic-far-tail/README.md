# Dynamic far-tail lemma draft

## Candidate statement

Let `Delta_j u` be periodic Littlewood-Paley blocks and let
`V_{>=m}` contain vortex-stretching triads whose largest and smallest
frequencies differ by more than `2^m`. For `m>=2`, the target estimate is

```text
|V_{>=m}(t)|
  <= C (2^(-m/2) + 2^(-3m/2)) Z(t)^(3/4) P(t)^(3/4),       (FT-1)
```

with `C` independent of the Galerkin cutoff. Here

```text
Z = ||grad u||_2^2,    P = ||Delta u||_2^2.
```

The C++ triad algebra and rational exponent checks are complete. Turning the
shell calculation below into a conventional cutoff-independent proof with all
projection constants stated is still open.

## Exact triad decomposition

Every separated triad has one unique low wave. The three roles are handled as
follows.

### Low advecting wave

For `k=p+q`, define

```text
T(p,q,k) = Re <u_k, i(q dot u_p)u_q>.
```

Reality and `p dot u_p=0` give the exact partner identity

```text
T(p,q,k)+T(p,-k,-q)=0,
|k|^2 T(p,q,k)+|q|^2 T(p,-k,-q)
    = (|k|^2-|q|^2)T(p,q,k).
```

Therefore

```text
abs(|k|^2-|q|^2) <= |p|(|k|+|q|),
```

which gains one low/high frequency factor.

### Low advected wave

When `q` is low, the multiplier `q dot u_p` already places the derivative on
the low wave. Relative to the cubic high-frequency weight this also gains one
low/high factor.

### Low target wave

When `k` is low, the enstrophy pairing carries the low weight `|k|^2`. This
gains two low/high factors and produces the stronger `2^(-3m/2)` shell decay.

`TriadTailEnvelope` verifies that these roles reconstruct the complete
nonlocal signed ledger and checks the coefficient-free Cauchy and frequency
bounds term by term.

## Shell summation

Write

```text
a_j = ||Delta_j u||_2.
```

Periodic Bernstein and Cauchy-Schwarz give

```text
sum_{j<=l-m} 2^(5j/2) a_j
    <= C 2^((l-m)/2) P^(1/2),

sum_{j<=l-m} 2^(7j/2) a_j
    <= C 2^(3(l-m)/2) P^(1/2).
```

The common high-frequency moment satisfies

```text
sum_l 2^(5l/2) a_l^2 <= Z^(3/4) P^(1/4).
```

The first inequality handles the low-advecting and low-advected roles; the
second handles the low-target role. Combining them yields (FT-1), subject to
writing out the bounded shell-overlap and torus Bernstein constants.

## Why a fixed gap is insufficient

The worst gap decay in (FT-1) gives

```text
|V_{>=m}|^4/(Z P^3) <= C 2^(-2m) Z^2.
```

The gap series is summable, but the energy identity controls `integral Z dt`,
not `integral Z^2 dt`. Therefore a fixed `m` does not close L4-A. This is
recorded as failed candidate F004.

## Moving-gap closure

Apply Young directly to (FT-1):

```text
|V_{>=m}| <= (nu/4) P + C nu^(-3) 2^(-2m) Z^3.
```

Choose the split dynamically,

```text
m(t) = m0 + ceil(log2(max(1,Z(t)))).
```

Then

```text
2^(-2m(t)) Z(t)^3 <= 2^(-2m0) Z(t),
```

so the geometric far tail contributes an absorbable palinstrophy term plus a
linear Gronwall term. `MovingGapController` verifies this integer/logarithmic
selection without using a bound on future enstrophy.

This conditionally closes the dynamically selected far tail, not the complete
Navier-Stokes estimate. All gaps below `m(t)` remain in a local/transition
block whose width is `O(log Z)`. Controlling that block without assuming the
desired regularity is the next lemma.

## Completion checklist

- [x] exact low-advecting commutator identity;
- [x] exact three-role finite Fourier reconstruction;
- [x] termwise frequency and amplitude inequalities;
- [x] rational Young and moving-gap exponent certificate;
- [ ] conventional Littlewood-Paley proof of (FT-1), including constants;
- [ ] cutoff-independent treatment of the `O(log Z)` transition/local block;
- [ ] insertion into the full Galerkin enstrophy inequality and compactness
      argument.
