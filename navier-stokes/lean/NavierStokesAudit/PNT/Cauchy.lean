module

public import NavierStokesAudit.PNT.Schur
public import Mathlib.Algebra.Order.BigOperators.Ring.Finset
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

/-! PNT-5/6/7: the selected channel, Cauchy alignment factor and squared majorant.
These identities require positive denominators. Cauchy-Schwarz supplies the two
pairing bounds; no smallness of the resulting majorant is inferred. -/
@[expose] public noncomputable section
namespace NavierStokesAudit.Cauchy
open scoped BigOperators

def selected (S s t Z P : ℝ) : ℝ := 3 * |S| * |s * t| / (2 * Z ^ 2 * P ^ 3)
def alignment (s t Z P B T : ℝ) : ℝ := s ^ 2 * t ^ 2 / (Z * P * B * T)
def majorant (S Z P B T : ℝ) : ℝ :=
  3 * |S| * Real.sqrt (Z * B * P * T) / (2 * Z ^ 2 * P ^ 3)
def majorantSq (S Z P B T : ℝ) : ℝ := 9 * S ^ 2 * B * T / (4 * Z ^ 3 * P ^ 5)

theorem finite_pairing {ι : Type*} [Fintype ι] (u v : ι → ℝ) :
    (∑ i, u i * v i) ^ 2 ≤ (∑ i, u i ^ 2) * ∑ i, v i ^ 2 :=
  Finset.sum_mul_sq_le_sq_mul_sq Finset.univ u v

theorem alignment_le_one (s t Z P B T : ℝ)
    (hZ : 0 < Z) (hP : 0 < P) (hB : 0 < B) (hT : 0 < T)
    (hs : s ^ 2 ≤ Z * B) (ht : t ^ 2 ≤ P * T) :
    alignment s t Z P B T ≤ 1 := by
  unfold alignment
  apply (div_le_one (by positivity)).mpr
  have h := mul_le_mul hs ht (sq_nonneg t) (by positivity : 0 ≤ Z * B)
  nlinarith [h]

theorem factorization (S s t Z P B T : ℝ)
    (hZ : 0 < Z) (hP : 0 < P) (hB : 0 < B) (hT : 0 < T) :
    selected S s t Z P = majorant S Z P B T * Real.sqrt (alignment s t Z P B T) := by
  have he : Z * B * P * T * alignment s t Z P B T = (s * t) ^ 2 := by
    unfold alignment
    field_simp
  have hr : Real.sqrt (Z * B * P * T) * Real.sqrt (alignment s t Z P B T) = |s * t| := by
    rw [← Real.sqrt_mul (by positivity), he, Real.sqrt_sq_eq_abs]
  unfold selected majorant
  rw [← hr]
  ring

theorem majorant_square (S Z P B T : ℝ)
    (hZ : 0 < Z) (hP : 0 < P) (hB : 0 ≤ B) (hT : 0 ≤ T) :
    majorant S Z P B T ^ 2 = majorantSq S Z P B T := by
  unfold majorant majorantSq
  rw [div_pow, mul_pow, mul_pow, sq_abs, Real.sq_sqrt (by positivity)]
  field_simp
  ring

theorem selected_le_majorant (S s t Z P B T : ℝ)
    (hZ : 0 < Z) (hP : 0 < P) (hB : 0 < B) (hT : 0 < T)
    (hs : s ^ 2 ≤ Z * B) (ht : t ^ 2 ≤ P * T) :
    selected S s t Z P ≤ majorant S Z P B T := by
  rw [factorization S s t Z P B T hZ hP hB hT]
  have ha := alignment_le_one s t Z P B T hZ hP hB hT hs ht
  have hsqrt : Real.sqrt (alignment s t Z P B T) ≤ 1 := by
    have h := Real.sqrt_le_sqrt ha
    simpa using h
  exact mul_le_of_le_one_right (by unfold majorant; positivity) hsqrt

end NavierStokesAudit.Cauchy
end
