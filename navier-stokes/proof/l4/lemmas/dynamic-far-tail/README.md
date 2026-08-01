# Dynamic far-tail lemma

## Cutoff-independent statement

Let `Delta_j u` be periodic Littlewood-Paley blocks and let
`V_{>=m}` contain vortex-stretching triads whose largest and smallest
frequencies differ by more than `2^m`. For `m>=2`, the hard-shell estimate is

```text
|V_{>=m}(t)|
  <= (C1 2^(-m/2) + C3 2^(-3m/2))
       Z(t)^(3/4) P(t)^(3/4),                              (FT-1)

C1 = 192 (2+2^(5/2)) sqrt(2),
C3 =  64 (2+2^(5/2)) sqrt(8/7).
```

with `C` independent of the Galerkin cutoff. Here

```text
Z = ||grad u||_2^2,    P = ||Delta u||_2^2.
```

The constants are independent of the Galerkin cutoff. The proof below uses
hard Fourier shells on the `2*pi` torus with normalized Haar measure. Because
`Delta u` is divergence-free and lies in the Galerkin space, the orthogonal
Leray and Galerkin projections disappear from its nonlinear pairing; they add
no operator constant.

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

## Periodic Fourier and shell constants

For the hard shell

```text
S_j = {k in Z^3 : 2^j <= |k| < 2^(j+1)},   j>=0,
```

the containing lattice cube has at most `64*2^(3j)` points. Parseval and
Cauchy-Schwarz therefore give

```text
sum_{k in S_j} |u_hat(k)| <= 8*2^(3j/2) a_j.
```

If a triad has gap at least two, its unique low magnitude is `L<H/4`. The
other high magnitude is greater than `3H/4`; hence the two high waves lie in
equal or adjacent hard shells.

Using the upper shell endpoints, the low-advecting commutator, low-advected,
and low-target frequency multipliers contribute respectively

```text
128 * 2^(5j/2) * 2^(2l),
 64 * 2^(5j/2) * 2^(2l),
 64 * 2^(7j/2) * 2^l
```

times the low and two high shell amplitudes. These constants include the
Fourier `l1` bound above. For nonnegative shell amplitudes and `s=5/2`,

```text
sum_{|r-t|<=1} 2^(s max(r,t)) a_r a_t
  <= (2+2^s) sum_l 2^(s l) a_l^2.
```

This follows from `2 a_r a_(r+1) <= a_r^2+a_(r+1)^2`; the coefficient of an
interior weighted square is exactly bounded by `2+2^s`.

## Shell summation

Write

```text
a_j = ||Delta_j u||_2.
```

After the torus Bernstein step, Cauchy-Schwarz and the exact geometric sums
give the cutoff-independent sequence inequalities

```text
sum_{j<=l-m} 2^(5j/2) a_j
    <= sqrt(2) 2^((l-m)/2) P^(1/2),

sum_{j<=l-m} 2^(7j/2) a_j
    <= sqrt(8/7) 2^(3(l-m)/2) P^(1/2).
```

The common high-frequency moment satisfies

```text
sum_l 2^(5l/2) a_l^2 <= Z^(3/4) P^(1/4).
```

The moment inequality has constant one by Holder interpolation. The first
low-shell inequality handles the low-advecting and low-advected roles; the
second handles the low-target role. Combining the three block constants, the
adjacent-shell constant, and these geometric sums gives exactly `C1` and `C3`
in (FT-1).

`PeriodicShellGeometry` checks every constant and the universal lattice cube
bound. `DyadicShellBounds` checks the scalar forms, and `PeriodicTailBound`
compares the complete signed Fourier tail against the final explicit right
side. These checks are regressions for the conventional calculation above;
the proof is the displayed finite Fourier and sequence inequalities.

## Why a fixed gap is insufficient

The worst gap decay in (FT-1) gives

```text
|V_{>=m}|^4/(Z P^3) <= C 2^(-2m) Z^2.
```

The gap series is summable, but the energy identity controls `integral Z dt`,
not `integral Z^2 dt`. Therefore a fixed `m` does not close L4-A. This is
recorded as failed candidate F004.

## Moving-gap closure

Set `A_m=C1*2^(-m/2)+C3*2^(-3m/2)`. Young's inequality with
conjugates `4/3` and `4` gives

```text
|V_{>=m}| <= (nu/4) P + 27/(4 nu^3) A_m^4 Z^3.
```

Choose the split dynamically,

```text
m(t) = m0 + ceil(log2(max(1,Z(t)))).
```

Using `(a+b)^4 <= 8(a^4+b^4)`, the moving-gap remainder satisfies

```text
27/(4 nu^3) A_m(t)^4 Z(t)^3 <= K(nu,m0) Z(t),

K(nu,m0) = 54 nu^(-3)
  * (C1^4 2^(-2m0) + C3^4 2^(-6m0)).
```

so the geometric far tail contributes an absorbable palinstrophy term plus a
linear Gronwall term. `MovingGapController` verifies the integer/logarithmic
selection without using a bound on future enstrophy, and `FarTailClosure`
checks the complete explicit Young remainder over logarithmically distributed
values of `Z`.

This closes the dynamically selected far-tail contribution uniformly in the
Galerkin cutoff. It does not close the complete Navier-Stokes estimate. All
gaps below `m(t)` remain in a local/transition block whose width is
`O(log Z)`. Controlling that block without assuming the desired regularity is
the next lemma.

## Completion checklist

- [x] exact low-advecting commutator identity;
- [x] exact three-role finite Fourier reconstruction;
- [x] termwise frequency and amplitude inequalities;
- [x] cutoff-independent scalar shell sums with explicit constants;
- [x] periodic Fourier mode count and adjacent high-shell estimate;
- [x] explicit cutoff-independent FT-1 constants;
- [x] rational Young and moving-gap exponent certificate;
- [x] explicit moving-gap linear-enstrophy remainder;
- [ ] cutoff-independent treatment of the `O(log Z)` transition/local block;
- [ ] insertion into the full Galerkin enstrophy inequality and compactness
      argument.
