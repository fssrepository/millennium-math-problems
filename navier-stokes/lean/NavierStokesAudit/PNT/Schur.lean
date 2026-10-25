module

public import Mathlib.Analysis.Real.Sqrt
public import Mathlib.Algebra.Order.BigOperators.Group.Finset
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

/-! Finite Gram-Schur estimate. Zero diagonal rows are retained with weight zero;
this is equivalent to omitting them, provided their Gram entries vanish. -/
@[expose] public noncomputable section
namespace NavierStokesAudit.Schur
open scoped BigOperators

variable {ι : Type*} [Fintype ι]

def normalized (g : ι → ι → ℝ) (d : ι → ℝ) (i j : ι) : ℝ :=
  |g i j| / Real.sqrt (d i * d j)

def row (g : ι → ι → ℝ) (d : ι → ℝ) (i : ι) : ℝ :=
  ∑ j, normalized g d i j

omit [Fintype ι] in
theorem normalized_nonneg (g : ι → ι → ℝ) (d : ι → ℝ) (i j : ι) :
    0 ≤ normalized g d i j := by
  exact div_nonneg (abs_nonneg _) (Real.sqrt_nonneg _)

theorem zero_row (g : ι → ι → ℝ) (d : ι → ℝ) (i : ι) (hi : d i = 0) :
    row g d i = 0 := by
  simp [row, normalized, hi]

omit [Fintype ι] in
theorem entry_bound (g : ι → ι → ℝ) (d : ι → ℝ)
    (hd : ∀ i, 0 ≤ d i)
    (hz : ∀ i j, d i = 0 ∨ d j = 0 → g i j = 0) (i j : ι) :
    g i j ≤ normalized g d i j * (d i + d j) / 2 := by
  by_cases hi : d i = 0
  · simp [normalized, hi, hz i j (Or.inl hi)]
  by_cases hj : d j = 0
  · simp [normalized, hj, hz i j (Or.inr hj)]
  have hip : 0 < Real.sqrt (d i) := Real.sqrt_pos.mpr (lt_of_le_of_ne (hd i) (Ne.symm hi))
  have hjp : 0 < Real.sqrt (d j) := Real.sqrt_pos.mpr (lt_of_le_of_ne (hd j) (Ne.symm hj))
  have he : |g i j| = normalized g d i j * (Real.sqrt (d i) * Real.sqrt (d j)) := by
    unfold normalized
    rw [Real.sqrt_mul (hd i)]
    field_simp
  have hm : Real.sqrt (d i) * Real.sqrt (d j) ≤ (d i + d j) / 2 := by
    nlinarith [sq_nonneg (Real.sqrt (d i) - Real.sqrt (d j)),
      Real.sq_sqrt (hd i), Real.sq_sqrt (hd j)]
  calc
    g i j ≤ |g i j| := le_abs_self _
    _ = normalized g d i j * (Real.sqrt (d i) * Real.sqrt (d j)) := he
    _ ≤ normalized g d i j * ((d i + d j) / 2) :=
      mul_le_mul_of_nonneg_left hm (normalized_nonneg g d i j)
    _ = _ := by ring

/-- PNT-11. R can be the maximum of the normalized absolute Gram rows.
For a Gram matrix, sum g is the squared norm of the sum of shell vectors. -/
theorem finite_bound (g : ι → ι → ℝ) (d : ι → ℝ) (R : ℝ)
    (hd : ∀ i, 0 ≤ d i)
    (hs : ∀ i j, g i j = g j i)
    (hz : ∀ i j, d i = 0 ∨ d j = 0 → g i j = 0)
    (hr : ∀ i, row g d i ≤ R) :
    (∑ i, ∑ j, g i j) ≤ R * ∑ i, d i := by
  have hw : ∀ i j, normalized g d i j = normalized g d j i := by
    intro i j
    simp only [normalized, hs i j, mul_comm]
  have sum_div (f : ι → ℝ) : (∑ i, f i / 2) = (∑ i, f i) / 2 := by
    simp only [div_eq_mul_inv, Finset.sum_mul]
  have swap : (∑ i, ∑ j, normalized g d i j * d j) =
      ∑ i, ∑ j, normalized g d i j * d i := by
    rw [Finset.sum_comm]
    congr 1
    funext i
    congr 1
    funext j
    rw [hw j i]
  calc
    (∑ i, ∑ j, g i j) ≤ ∑ i, ∑ j, normalized g d i j * (d i + d j) / 2 := by
      apply Finset.sum_le_sum; intro i _
      apply Finset.sum_le_sum; intro j _
      exact entry_bound g d hd hz i j
    _ = ∑ i, row g d i * d i := by
      simp only [mul_add, add_div, Finset.sum_add_distrib, sum_div]
      rw [swap]
      simp only [row, ← Finset.sum_mul]
      ring
    _ ≤ ∑ i, R * d i := by
      apply Finset.sum_le_sum; intro i _
      exact mul_le_mul_of_nonneg_right (hr i) (hd i)
    _ = R * ∑ i, d i := by rw [Finset.mul_sum]

variable {κ : Type*} [Fintype κ]

def vectorGram (v : ι → κ → ℝ) (i j : ι) : ℝ := ∑ k, v i k * v j k
def vectorDiagonal (v : ι → κ → ℝ) (i : ι) : ℝ := ∑ k, v i k ^ 2

omit [Fintype ι] in
theorem vector_zero (v : ι → κ → ℝ) (i : ι) (hi : vectorDiagonal v i = 0) :
    ∀ k, v i k = 0 := by
  intro k
  have hz : ∀ k ∈ Finset.univ, v i k ^ 2 = 0 :=
    (Finset.sum_eq_zero_iff_of_nonneg (fun k _ => sq_nonneg (v i k))).mp hi
  have h := hz k (Finset.mem_univ k)
  nlinarith [sq_nonneg (v i k)]

/-- Gram specialization: the matrix hypotheses follow from the vectors,
including every zero-diagonal row. This is the finite shell-vector inequality. -/
theorem vector_bound (v : ι → κ → ℝ) (R : ℝ)
    (hr : ∀ i, row (vectorGram v) (vectorDiagonal v) i ≤ R) :
    (∑ k, (∑ i, v i k) ^ 2) ≤ R * ∑ i, ∑ k, v i k ^ 2 := by
  have hd : ∀ i, 0 ≤ vectorDiagonal v i := by
    intro i; exact Finset.sum_nonneg (fun k _ => sq_nonneg (v i k))
  have hs : ∀ i j, vectorGram v i j = vectorGram v j i := by
    intro i j; simp only [vectorGram, mul_comm]
  have hz : ∀ i j, vectorDiagonal v i = 0 ∨ vectorDiagonal v j = 0 →
      vectorGram v i j = 0 := by
    intro i j h
    rcases h with hi | hj
    · simp [vectorGram, vector_zero v i hi]
    · simp [vectorGram, vector_zero v j hj]
  have he : (∑ k, (∑ i, v i k) ^ 2) = ∑ i, ∑ j, vectorGram v i j := by
    simp only [sq, Finset.sum_mul, Finset.mul_sum, vectorGram]
    rw [Finset.sum_comm]
    apply Finset.sum_congr rfl
    intro i _
    rw [Finset.sum_comm]
    simp only [mul_comm]
  rw [he]
  exact finite_bound (vectorGram v) (vectorDiagonal v) R hd hs hz hr

end NavierStokesAudit.Schur
end
