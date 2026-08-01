# Local triad energy symmetry

This directory records an exact partial lemma for the remaining gap-zero
block. It is independent of the Galerkin cutoff and does not use numerical
boundedness.

## LS-1

Fix one unordered Fourier triad `a+b+c=0`. Let `T_a`, `T_b`, and `T_c` be its
complete nonlinear kinetic-energy transfers into the three target modes,
including both orders of the two input modes. Incompressibility and the
skew-symmetry of advection give the triadwise identity

```text
T_a + T_b + T_c = 0.
```

Its contribution to enstrophy stretching is

```text
V_triad = |a|^2 T_a + |b|^2 T_b + |c|^2 T_c.
```

For any real reference `r`, the energy identity gives

```text
V_triad = (|a|^2-r)T_a + (|b|^2-r)T_b + (|c|^2-r)T_c.
```

Choosing `r=min(|a|^2,|b|^2,|c|^2)` proves

```text
|V_triad|
 <= (max(|a|^2,|b|^2,|c|^2)-min(|a|^2,|b|^2,|c|^2))
    (|T_a|+|T_b|+|T_c|).                                  (LS-1)
```

In particular, an equal-length triad has exactly zero enstrophy transfer even
when its individual ordered interactions are nonzero. `LS-1` is a conventional
finite algebraic identity and remains valid after summing any set of complete
triads.

## Machine verification

`LocalTriadSymmetrizer` groups the ordered convolution entries by the exact
signed key `sort(-k,p,q)`, reconstructs the gap-zero enstrophy ledger, and
checks both the triad energy identity and `LS-1`. The deterministic self-test
currently reports

```text
maximum relative energy-cancellation residual = 8.70e-19
local-ledger reconstruction residual           = 1.16e-20
maximum LS-1 envelope ratio                     = 4.99e-1
```

The K3 heterochiral endpoint gives the following spread decomposition:

```text
relative squared-frequency spread   signed V       absolute grouped V
equal                               -3.26e-21       6.66e-21
(0,1/4]                              1.03e-6        1.03e-6
(1/4,1/2]                           -2.07019e-1     2.07019e-1
(1/2,3/4]                           -1.99194e-2     1.99956e-2
```

The replayable detailed certificate is
`../../analysis/helical-heterochiral-cutoff/K3.json`.

## Broad-spread adversary

The exact sector objective and its checkpointed RK4 adjoint accept the same
frequency-spread mask. Here

```text
narrow: 4(max |k_i|^2-min |k_i|^2) <= max |k_i|^2,
broad:  the strict complementary inequality inside gap zero.
```

The broad masked gradient agrees with central differences at approximately
`1.1e-11`. Eight steps on each of twelve parallel K3 starts increase its
trajectory objective from `4.39053937751e-6` to `4.39098016662e-6`. On the
winning state, evaluation without further optimization gives

```text
all local heterochiral     4.39089290346e-6
broad-spread               4.39098016662e-6
narrow-spread              2.71299915671e-27
equal-frequency            8.00967140430e-89
```

The fourth-power objectives are not additive, but the scale separation is
decisive: the current local extremum is a broad-spread extremum. Its same-state
cutoff scan is

```text
K:       3                  4                  5                  6
J: 4.39098016662e-6   4.39098438285e-6   4.39098438272e-6   4.39098438272e-6
```

with final adjacent relative difference `2.66e-15`. See
`../../adversary/helical-heterochiral-broad-spread-cutoff-K3-K6.json` and
`../../analysis/helical-heterochiral-broad-spread/K3.json`.

## What remains open

`LS-1` removes every equal-frequency local triad, but it does not yet close
L4. The dominant adversarial transfer lies in the `(1/4,1/2]` spread bin, where
the frequency difference is still comparable to the squared frequency. Thus
the identity supplies no negative power of the active scale on the entire
local block. A successful continuation must combine `LS-1` with a dynamic
bound on broad-spread local triads, or find a different time-integrated
mechanism that supplies the missing `Z^(-1/2)` depletion.
