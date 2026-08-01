# Projected normalization residual

This identity combines three of the five local quartet entries before any
inequality is applied.

Write

```text
S=<A u,B_L>,  Z=<u,A u>,  P=<A u,A u>,  h=-B_L,
B_L=B_local(u,u).
```

The outer Frechet slot and the `Z'` and `P'` chain-rule terms in
`C=S^4/(ZP^3)` have the common coefficient `4S^3/(ZP^3)`. Using

```text
Z'[h]=2<A u,h>,   P'[h]=2<A u,A h>,
```

their exact combination is

```text
4S^3/(ZP^3) <A h,R>,

R = B_L - S u/(2Z) - 3S A u/(2P).                    (LQR-1)
```

For `h=-B_L`, self-adjointness of `A` and `<A u,B_L>=S` expand the raw
pairing as

```text
<A h,R>
 = -||A^(1/2)B_L||_2^2
   + S^2/(2Z)
   + 3S <A B_L,A u>/(2P).                            (LQR-2)
```

Completing the square in the last two `A^(1/2)B_L` terms gives the equivalent
identity

```text
beta = 3S/(4P),

<A h,R>
 = -||A^(1/2)(B_L-beta A u)||_2^2
   + S^2/(2Z)
   + beta^2 H3,

H3=||A^(3/2)u||_2^2.                                 (LQR-3)
```

After LQR-1, the complete normalized local nonlinear contribution has only
three entries:

```text
normalized projected pairing + normalized advecting slot
                             + normalized advected slot.
```

On the optimized K6 state they are

```text
projected pairing             +1.9943549392e-3
advecting slot                -6.1297970097e-5
advected slot                 -1.2458880655e-3
total                         +6.8716890362e-4.
```

The two independent forms of the raw projected pairing are

```text
direct                         -0.521682741616
expanded LQR-2                 -0.521682741616
completed LQR-3                -0.521682741616
```

The expansion, square-completion, and normalized full-quartet residuals are
`1.56e-18`, `5.20e-19`, and `3.85e-18`, respectively.

The negative square in LQR-3 is not automatically favorable: its multiplier
contains `S^3` and changes sign when the signed local stretching changes sign.
Thus LQR-3 is an exact reduction, not a sign proof. A successful estimate must
be sign-aware and must combine the projected pairing with at least the two
remaining transport variations; separately taking their absolute values
would restore an uncontrolled `H3` contribution.
