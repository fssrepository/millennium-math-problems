# Projective self/cross quartet decomposition

This note isolates the genuinely open part of the projective remainder. It
does not prove a cutoff-independent cross-ray bound or the full shifted local
density lemma.

## Exact decomposition

Write the selected local advection as

```text
B(u,u) = sum_sigma B_sigma(u,u) = sum_sigma b_sigma,
c_sigma = B_sigma(u,Au),
s_sigma = <Au,b_sigma>,
t_sigma = <A b_sigma,Au>.
```

For one primitive projective ray define its self-quartet

```text
J_sigma = -<b_sigma,A b_sigma>
          +<b_sigma,c_sigma>
          -<Au,B_sigma(b_sigma,u)>
          +s_sigma^2/(2Z)+3s_sigma t_sigma/(2P).     (PCQ-1)
```

Expanding every occurrence of `B=sum_sigma B_sigma` in the exact selected
bracket gives

```text
J_selected = sum_sigma J_sigma + J_cross.           (PCQ-2)
```

Here `J_cross` contains exactly the unequal-ray pairs:

```text
J_cross = -sum_(sigma!=tau)<b_sigma,A b_tau>
          +sum_(sigma!=tau)<b_sigma,c_tau>
          -sum_(sigma!=tau)<Au,B_tau(b_sigma,u)>
          +[(sum s_sigma)^2-sum s_sigma^2]/(2Z)
          +3[(sum s_sigma)(sum t_sigma)
             -sum s_sigma t_sigma]/(2P).            (PCQ-3)
```

`ProjectiveQuarticDiagonalKernel` evaluates `sum J_sigma` directly. Its
sparse reverse pass differentiates PCQ-1 without allocating a full Fourier
field for every ray. The central-difference error of the resulting cross-only
power-one objective is `2.21e-10` on the self-test state.

## What the decomposition changes

On the original power-one maximizers, the self part dominates:

```text
K       full power one       self rays          cross rays
1       0.00600652           0.00601385         -0.00000733
2       0.00629048           0.00616128          0.00012920
3       0.00634032           0.00616229          0.00017803
4       0.00635502           0.00614804          0.00020698
5       0.00636095           0.00615762          0.00020333
6       0.00636381           0.00614329          0.00022051
7       0.00636477           0.00614715          0.00021762
```

At K7 the effective self-ray count is `1.076`, and the dominant ray carries
`96.36%` of the absolute self sum. After removing the proved double/triple
families and the exact `(1,2,3)` signature, the K5 tail still decomposes as

```text
full=0.00256213,  self=0.00202885,  cross=0.00053328.
```

Thus the largest known power-one states are primarily fixed-ray states. This
does not by itself sum the self bounds uniformly over all primitive shapes.

## Cross-only adversary

The exact-gradient objective

```text
|J_cross S_full|^2/(Z^4 P^4)                        (PCQ-4)
```

subtracts every self-quartet before optimization. Twelve-start searches give
the following square roots of PCQ-4:

```text
K                 1          2          3          4          5
|J_cross S|/Z^2P^2 0.0004531  0.0010594  0.0011000  0.0011144  0.0011166
```

A small-step K6 continuation reaches `0.00111672`; exact evaluate-only lifts
through K8 retain the same value. On the K5 cross winner, the self and cross
pieces have opposite signs:

```text
full=+0.00059465,  self=+0.00171121,  cross=-0.00111656.
```

The flat branch and the self/cross cancellation are finite evidence. They do
not exclude another high-frequency cross branch. The remaining analytic task
is a cutoff-independent estimate for PCQ-3 after multiplication by
`S_full/(Z^2P^2)`, together with a uniform summation of PCQ-1.

## Cross attribution and fixed-core closure

The exact per-ray attribution ledger writes each ray's cross contribution as
its full primary-ray attribution minus its self quartet. This is an
attribution, not a pairwise ray-by-ray matrix. On the K7 power-one winner the
cross block has `2.103` effective rays and the leading `(1,2,3)` ray carries
`64.96%`. On the K6 cross-only winner it has `3.583` effective rays and the
leading ray carries `45.78%`; the first four attributions are `(1,2,3)`,
`(2,3,5)`, `(1,1,1)`, and `(1,3,4)`.

This concentration leads to an analytic reduction stronger than a finite
diagnostic. For any fixed finite ray set `F`, its union incidence is still
`O_F(K)`. Applying the fixed-ray bilinear argument to the unexpanded operator
`B_F=sum_(sigma in F)B_sigma` closes the complete internal self+cross quartet
with a cutoff-independent constant and frequency gain `-1/2`. See
[FINITE_PROJECTIVE_FAMILY.md](FINITE_PROJECTIVE_FAMILY.md).

Therefore PCQ-3 is no longer open inside any fixed finite core. The remaining
statement is precisely the coupling between such a core and the growing
projective tail, plus the tail's internal block. The constant is not uniform
when the core itself grows with the cutoff.

## Reproduction

```bash
./build/navier_stokes_lab local-sld-projective-quartic-cross \
  --state proof/l4/states/local-shape-power/remainder-p1-K7/K7.tsv \
  --threads 12 --top 64 \
  --certificate proof/l4/analysis/shifted-local-density/remainder-quartet/K7-power-one-quartic-cross.json

./build/navier_stokes_lab local-closure-adversary \
  --objective projective-cross-power-ratio \
  --selection double-triple-remainder \
  --min-cutoff 5 --max-cutoff 5 \
  --restarts 12 --workers 12 --iterations 8 \
  --method lbfgs --backend direct \
  --warm-state proof/l4/states/local-projective-cross-power/double-triple-remainder-K1-K5/restarts/K4/R000.tsv \
  --certificate proof/l4/adversary/shifted-local-density/projective-cross-power/K5-warm.json \
  --state-dir proof/l4/states/local-projective-cross-power/K5-warm

./build/navier_stokes_lab local-closure-adversary \
  --objective projective-cross-power-ratio \
  --selection double-triple-remainder \
  --min-cutoff 8 --max-cutoff 8 --restarts 1 --workers 12 \
  --iterations 0 --backend direct \
  --warm-state proof/l4/states/local-projective-cross-power/K7-lift-evaluate/K7.tsv \
  --certificate proof/l4/adversary/shifted-local-density/projective-cross-power/K8-lift-evaluate.json \
  --state-dir proof/l4/states/local-projective-cross-power/K8-lift-evaluate
```

Artifacts:

- [`../../analysis/shifted-local-density/remainder-quartet/K7-power-one-quartic-cross.json`](../../analysis/shifted-local-density/remainder-quartet/K7-power-one-quartic-cross.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K5-power-one-tail-quartic-cross.json`](../../analysis/shifted-local-density/remainder-quartet/K5-power-one-tail-quartic-cross.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K5-cross-power-winner-quartic-cross.json`](../../analysis/shifted-local-density/remainder-quartet/K5-cross-power-winner-quartic-cross.json)
- [`../../adversary/shifted-local-density/projective-cross-power/K5-warm.json`](../../adversary/shifted-local-density/projective-cross-power/K5-warm.json)
- [`../../adversary/shifted-local-density/projective-cross-power/K6-small-step.json`](../../adversary/shifted-local-density/projective-cross-power/K6-small-step.json)
- [`../../adversary/shifted-local-density/projective-cross-power/K8-lift-evaluate.json`](../../adversary/shifted-local-density/projective-cross-power/K8-lift-evaluate.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K7-power-one-cross-attribution.json`](../../analysis/shifted-local-density/remainder-quartet/K7-power-one-cross-attribution.json)
- [`../../analysis/shifted-local-density/remainder-quartet/K6-cross-power-winner-cross-attribution.json`](../../analysis/shifted-local-density/remainder-quartet/K6-cross-power-winner-cross-attribution.json)
