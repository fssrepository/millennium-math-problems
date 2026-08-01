# Local squared-length signature families

This directory proves a cutoff-uniform estimate for every fixed local
squared-length signature and isolates one precise remaining condition for the
broad local block. It extends the equal-length orthogonal result from
`(r,r,2r)` to arbitrary ordered triples `(r,s,t)`.

## Fixed-signature lattice degree

Fix a nonzero input wave `p` and an ordered signature

```text
|p|^2 = r,       |q|^2 = s,       |p+q|^2 = t.
```

The last two equations imply the plane constraint

```text
2 p dot q = t-r-s.
```

Choose a nonzero component of `p` and eliminate the corresponding component
of `q`. Substitution into `|q|^2=s` leaves a positive-definite quadratic in
the other two integer coordinates. After one coordinate is fixed, there are
at most two values of the other. In a component cutoff-`K` cube,

```text
degree_p(r,s,t) <= 2(2K+1).                              (LSF-1)
```

Fixing the target and writing `q=k-p` gives the same plane/sphere
intersection and the same target-degree bound. Ordering or symmetrizing the
three lengths changes only a finite constant.

## Square-summed transfer lemma

Let `V_sigma` be the complete-triad enstrophy transfer for one local
signature in a shell of frequency `K`. The complete-triad cancellation LS-1
supplies two frequency powers and advection supplies one. For scalar Fourier
amplitudes `a`, `b`, and `c`, set

```text
X_sigma = sum_edges(sigma) a_p^2 b_q^2,
Y_sigma = sum_edges(sigma) c_k^2.
```

Cauchy--Schwarz and LSF-1 give

```text
|V_sigma|^2 <= C K^6 X_sigma Y_sigma,
Y_sigma <= C K ||c||_2^2.
```

Every ordered interaction belongs to exactly one ordered signature, hence

```text
sum_sigma X_sigma <= ||a||_2^2 ||b||_2^2.
```

After summing before taking the square root,

```text
(sum_sigma |V_sigma|^2)^(1/2)
    <= C K^(7/2) E_K^(3/2).                              (LSF-2)
```

This is stronger than applying a separate bound to every signature and then
counting all arithmetically possible signatures.

## Exact effective-count reduction

For a nonzero transfer vector define its inverse participation number

```text
N_eff = (sum_sigma |V_sigma|)^2 / sum_sigma |V_sigma|^2.
```

This is an identity, not a statistical model. Combining it with LSF-2 yields

```text
|sum_sigma V_sigma|
 <= C sqrt(N_eff) K^(7/2) E_K^(3/2).                    (LSF-3)
```

If a trajectory-uniform analytical estimate proves
`N_eff(K,t) <= C K^mu`, then the transfer frequency power is
`7/2+mu/2`. Relative to `K^4 E_K` viscous dissipation, the remaining factor is

```text
K^(-1/2+mu/2) E(0)^(1/2).
```

Consequently every `mu<1` closes at high frequency. The case `mu=1` is the
critical boundary, while the dense `mu=2` count restores the standard
supercritical `K^(9/2)` estimate.

There is a sharper exact reduction that does not discard inter-signature
signs. Define

```text
A_sig = |sum_sigma V_sigma| / (sum_sigma |V_sigma|^2)^(1/2).
```

Then LSF-2 directly gives

```text
|sum_sigma V_sigma| <= C A_sig K^(7/2) E_K^(3/2).        (LSF-4)
```

If `A_sig=O(K^alpha)`, viscosity absorbs the high shells whenever
`alpha<1/2`. This is strictly more informative than replacing the signed sum
by `sum |V_sigma|` and bounding `N_eff`.

## Machine certificates

`LocalSignatureGeometry` performs exact exponent arithmetic and enumerates
the lattice degree as a regression check. Through cutoff five it reports

```text
K   ordered signatures   max input degree   max target degree   2(2K+1)
1           10                   4                   4                6
2          117                   8                   8               10
3          659                   8                   8               14
4         2633                  12                  12               18
5         7730                  12                  12               22
```

The machine-readable artifact is
[`local-signature-certificate.json`](local-signature-certificate.json).
Reproduce it with

```bash
./build/navier_stokes_lab local-signature-certificate \
  --max-cutoff 5 \
  --certificate proof/l4/lemmas/local-signature-families/local-signature-certificate.json
```

`LocalTriadSymmetrizer` also computes `N_eff` directly from replayed states.
For the current K3 broad heterochiral endpoint it finds

```text
coherent signatures above numerical threshold: 106
N_eff:                                      1.1902068038
dominant fraction of sum |V_sigma|:         0.912689851446
signed signature amplification A_sig:       1.09095551858
```

These values are stored in
[`../../../analysis/helical-heterochiral-broad-spread/K3.json`](../../../analysis/helical-heterochiral-broad-spread/K3.json).
They are numerical evidence only; LSF-1--LSF-3 do not depend on them.

## Remaining lemma

The local broad-spread obstruction is now reduced to a concrete analytical
statement: prove a trajectory-integrated or viscosity-weighted sublinear
bound on the effective coherent signature count, or obtain an equivalent
cancellation that bounds the left side of LSF-3 directly. A 12-worker
adversarial search reaches `N_eff=1748.62` at K6 for arbitrary flat-spectrum
unit-energy states, so a pointwise sublinear estimate based only on energy and
incompressibility fails the computational acceptance test. The K3 broad
trajectory value near one remains a strong dynamic pattern, but it cannot be
promoted to an unrestricted algebraic lemma. See
[`../../analysis/local-signature-adversary/README.md`](../../analysis/local-signature-adversary/README.md).

The sharper LSF-4 candidate survives the same finite search: its K6 maxima are
`2.33549` on flat states and `2.73471` on outer-half-flat states, with fitted
exponents `0.48310` and `0.25607`. The next proof target is a
cutoff-independent `A_sig=O(K^alpha)` estimate with `alpha<1/2`, or its
time-integrated analogue.
