module

public import Mathlib.Basic.Real.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

/-!
Exact scalar formulas from the PNT normalization objective.
The statements in this module are algebraic implications. Fourier realizability
and a nonzero numerator are supplied by the separate concrete witness module.
-/

@[expose] public noncomputable section

namespace NavierStokesAudit

structure PNTData where
  energy : ℝ
  Z : ℝ
  P : ℝ
  S : ℝ
  sc : ℝ
  st : ℝ
  tc : ℝ
  tt : ℝ
  B : ℝ
  D : ℝ
  R : ℝ

def PNTData.product (d : PNTData) : ℝ := d.sc * d.tt + d.st * d.tc + d.st * d.tt

/-- The absolute open objective PNT-2, including its factor 3/2. -/
def PNTData.openValue (d : PNTData) : ℝ :=
  3 * |d.S| * |d.product| / (2 * d.Z ^ 2 * d.P ^ 3)

/-- PNT-7 when R = 1 and D = T2; the PNT-12 Schur formula otherwise. -/
def PNTData.schurValue (d : PNTData) : ℝ :=
  9 * d.S ^ 2 * d.B * d.R * d.D / (4 * d.Z ^ 3 * d.P ^ 5)

/-- Amplitude action on the scalar quantities, for a fixed spectral geometry. -/
def PNTData.scale (a : ℝ) (d : PNTData) : PNTData where
  energy := a ^ 2 * d.energy
  Z := a ^ 2 * d.Z
  P := a ^ 2 * d.P
  S := a ^ 3 * d.S
  sc := a ^ 3 * d.sc
  st := a ^ 3 * d.st
  tc := a ^ 3 * d.tc
  tt := a ^ 3 * d.tt
  B := a ^ 4 * d.B
  D := a ^ 4 * d.D
  R := d.R

theorem core_tail_factorization (sc st tc tt : ℝ) :
    (sc + st) * (tc + tt) - sc * tc = sc * tt + st * tc + st * tt := by
  ring

theorem core_tail_two_term (sc st tc tt : ℝ) :
    sc * tt + st * tc + st * tt = (sc + st) * tt + st * tc := by
  ring

theorem product_scale (d : PNTData) (a : ℝ) :
    (d.scale a).product = a ^ 6 * d.product := by
  simp only [PNTData.scale, PNTData.product]
  ring

theorem open_scale (d : PNTData) (a : ℝ) (ha : 0 < a)
    (hZ : d.Z ≠ 0) (hP : d.P ≠ 0) :
    (d.scale a).openValue = d.openValue / a := by
  rw [PNTData.openValue, product_scale]
  simp only [PNTData.scale, abs_mul, abs_of_pos (pow_pos ha 3),
    abs_of_pos (pow_pos ha 6), PNTData.openValue]
  field_simp

theorem schur_scale (d : PNTData) (a : ℝ) (ha : a ≠ 0)
    (hZ : d.Z ≠ 0) (hP : d.P ≠ 0) :
    (d.scale a).schurValue = d.schurValue / a ^ 2 := by
  simp only [PNTData.schurValue, PNTData.scale]
  field_simp

/-- Construct a positive amplitude which defeats a proposed bound. -/
theorem inverse_exceeds (q M : ℝ) (hq : 0 < q) :
    ∃ a : ℝ, 0 < a ∧ a < 1 ∧ M < q / a ∧ M < q / a ^ 2 := by
  let a := q / (q + |M| + 1)
  have hd : 0 < q + |M| + 1 := by positivity
  have ha : 0 < a := div_pos hq hd
  have ha1 : a < 1 := (div_lt_one hd).2 (by linarith [abs_nonneg M])
  have hv : q / a = q + |M| + 1 := by
    dsimp [a]
    field_simp
  have hm : M < q / a := by rw [hv]; linarith [le_abs_self M]
  have hsq : a ^ 2 ≤ a := by nlinarith
  have hle : q / a ≤ q / a ^ 2 :=
    div_le_div_of_nonneg_left (le_of_lt hq) (pow_pos ha 2) hsq
  exact ⟨a, ha, ha1, hm, lt_of_lt_of_le hm hle⟩

/-- Conditional scalar result: a positive base value excludes a uniform bound. -/
theorem open_unbounded (d : PNTData) (hZ : d.Z ≠ 0) (hP : d.P ≠ 0)
    (hq : 0 < d.openValue) :
    ¬ ∃ C : ℝ, ∀ a : ℝ, 0 < a → (d.scale a).openValue ≤ C := by
  rintro ⟨C, hC⟩
  obtain ⟨a, ha, _, hgt, _⟩ := inverse_exceeds d.openValue C hq
  have hbound := hC a ha
  rw [open_scale d a ha hZ hP] at hbound
  exact (not_lt_of_ge hbound) hgt

/-- The squared majorant has the same obstruction, with degree minus two. -/
theorem schur_unbounded (d : PNTData) (hZ : d.Z ≠ 0) (hP : d.P ≠ 0)
    (hq : 0 < d.schurValue) :
    ¬ ∃ C : ℝ, ∀ a : ℝ, 0 < a → (d.scale a).schurValue ≤ C := by
  rintro ⟨C, hC⟩
  obtain ⟨a, ha, _, _, hgt⟩ := inverse_exceeds d.schurValue C hq
  have hbound := hC a ha
  rw [schur_scale d a (ne_of_gt ha) hZ hP] at hbound
  exact (not_lt_of_ge hbound) hgt

theorem fixed_energy_iff (d : PNTData) (hE : 0 < d.energy)
    (a : ℝ) (ha : 0 < a) :
    (d.scale a).energy = d.energy ↔ a = 1 := by
  simp only [PNTData.scale]
  constructor
  · intro h
    have hs : a ^ 2 = 1 := by nlinarith
    nlinarith
  · rintro rfl
    ring

theorem amplitude_degrees :
    (3 + 3 + 3 - (2 * 2 + 3 * 2) : ℤ) = -1 ∧
    (2 * 3 + 4 + 0 + 4 - (3 * 2 + 5 * 2) : ℤ) = -2 := by
  norm_num

end NavierStokesAudit

end
