module

public import NavierStokesAudit.Fourier.Cube
public import NavierStokesAudit.PNT.Schur
import Mathlib.Algebra.BigOperators.Fin
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

/-! PNT-7/8/12/14 at H=1 on the complete cutoff-one Fourier cube.
All nonzero interactions have primitive height 3, so shell 2 is the only
possible nonzero tail shell. The zero rows 0 and 1 are included explicitly;
CubeWitness.other_shell_advection proves that every later shell is zero too. -/
@[expose] public noncomputable section
namespace NavierStokesAudit.Majorants
open FourierWitness
open scoped BigOperators

abbrev b (u : State) := CubeWitness.advection CubeWitness.selected u
abbrev shell (u : State) (j : Fin 3) :=
  CubeWitness.advection (CubeWitness.shellSelection j.val) u

def shellGram (u : State) (i j : Fin 3) : ℝ := CubeWitness.gram (shell u i) (shell u j)
def diagonal (u : State) (j : Fin 3) : ℝ := CubeWitness.h2 (shell u j)
def diagonalSum (u : State) : ℝ := ∑ j, diagonal u j
def gramRow (u : State) (i : Fin 3) : ℝ := Schur.row (shellGram u) (diagonal u) i
def maxRow (u : State) : ℝ := max (gramRow u 0) (max (gramRow u 1) (gramRow u 2))

def data (u : State) : PNTData :=
  { FourierWitness.data u with
    B := CubeWitness.h1 (b u)
    D := diagonalSum u
    R := maxRow u }

/-- PNT-7, before the height factor in PNT-8. -/
def cauchy (u : State) : ℝ :=
  9 * (FourierWitness.data u).S ^ 2 * CubeWitness.h1 (b u) * CubeWitness.h2 (b u) /
    (4 * enstrophy u ^ 3 * palinstrophy u ^ 5)

/-- PNT-14 at a specified Gram row and height threshold one. -/
def fixedRow (u : State) (i : Fin 3) : ℝ :=
  9 * (FourierWitness.data u).S ^ 2 * CubeWitness.h1 (b u) * gramRow u i * diagonalSum u /
    (4 * enstrophy u ^ 3 * palinstrophy u ^ 5)

def compensated (H q : ℝ) : ℝ := Real.sqrt H * q

theorem shell_zero (u : State) : shell u 0 = fun _ _ => 0 :=
  CubeWitness.other_shell_advection 0 (by decide) u

theorem shell_one (u : State) : shell u 1 = fun _ _ => 0 :=
  CubeWitness.other_shell_advection 1 (by decide) u

theorem shell_two (u : State) : shell u 2 = b u := by
  change CubeWitness.advection (CubeWitness.shellSelection 2) u = b u
  rw [CubeWitness.shell_two]

theorem tail_sum (u : State) : (fun k j => ∑ i : Fin 3, shell u i k j) = b u := by
  funext k j
  simp [Fin.sum_univ_succ, shell_zero, shell_one, shell_two]

theorem gram_velocity (a : ℝ) (i j : Fin 3) :
    shellGram (velocity a) i j = if i = 2 ∧ j = 2 then 50 * a ^ 4 else 0 := by
  have hh : CubeWitness.gram (b (velocity a)) (b (velocity a)) = 50 * a ^ 4 :=
    CubeWitness.h2_velocity a
  have hz (c : CubeWitness.Output) : CubeWitness.gram (fun _ _ => 0) c = 0 := by
    simp [CubeWitness.gram, CubeWitness.sum27, dot]
  have hz' (c : CubeWitness.Output) : CubeWitness.gram c (fun _ _ => 0) = 0 := by
    simp [CubeWitness.gram, CubeWitness.sum27, dot]
  fin_cases i <;> fin_cases j <;>
    simp [shellGram, shell_zero, shell_one, shell_two, hh, hz, hz']

theorem diagonal_velocity (a : ℝ) (j : Fin 3) :
    diagonal (velocity a) j = if j = 2 then 50 * a ^ 4 else 0 := by
  change shellGram (velocity a) j j = _
  rw [gram_velocity]
  simp

theorem diagonal_sum_velocity (a : ℝ) : diagonalSum (velocity a) = 50 * a ^ 4 := by
  simp [diagonalSum, diagonal_velocity]

theorem gram_row_velocity (a : ℝ) (ha : 0 < a) (i : Fin 3) :
    gramRow (velocity a) i = if i = 2 then 1 else 0 := by
  have hd : 0 < 50 * a ^ 4 := by positivity
  fin_cases i <;>
    simp [gramRow, Schur.row, Schur.normalized, gram_velocity, diagonal_velocity,
      Fin.sum_univ_succ, Real.sqrt_mul_self (le_of_lt hd), abs_of_pos hd, ne_of_gt hd]

theorem max_row_velocity (a : ℝ) (ha : 0 < a) : maxRow (velocity a) = 1 := by
  norm_num [maxRow, gram_row_velocity a ha]

def baseData : PNTData := ⟨8, 18, 46, 2, 0, 2, 0, 6, 22, 50, 1⟩

theorem data_velocity (a : ℝ) (ha : 0 < a) : data (velocity a) = baseData.scale a := by
  unfold data
  rw [FourierWitness.data_velocity, diagonal_sum_velocity, max_row_velocity a ha]
  simp [b, CubeWitness.h1_velocity, FourierWitness.baseData, baseData, PNTData.scale, mul_comm]

theorem base_value : baseData.schurValue = 275 / 33366002112 := by
  norm_num [baseData, PNTData.schurValue]

theorem schur_value (a : ℝ) (ha : 0 < a) :
    (data (velocity a)).schurValue = 275 / (33366002112 * a ^ 2) := by
  rw [data_velocity a ha, schur_scale baseData a (ne_of_gt ha)
    (by norm_num [baseData]) (by norm_num [baseData]), base_value]
  ring

theorem cauchy_value (a : ℝ) (ha : 0 < a) :
    cauchy (velocity a) = 275 / (33366002112 * a ^ 2) := by
  simp only [cauchy, FourierWitness.data_velocity, FourierWitness.baseData, PNTData.scale,
    b, CubeWitness.h1_velocity, CubeWitness.h2_velocity, enstrophy_velocity, palinstrophy_velocity]
  have hn := ne_of_gt ha
  field_simp
  ring

theorem fixed_row_value (a : ℝ) (ha : 0 < a) :
    fixedRow (velocity a) 2 = 275 / (33366002112 * a ^ 2) := by
  simp only [fixedRow, FourierWitness.data_velocity, FourierWitness.baseData, PNTData.scale,
    b, CubeWitness.h1_velocity, diagonal_sum_velocity, gram_row_velocity a ha,
    enstrophy_velocity, palinstrophy_velocity, ↓reduceIte]
  have hn := ne_of_gt ha
  field_simp
  ring

/-- H=1 is an allowed height threshold; its compensation factor is exactly one. -/
theorem height_one (q : ℝ) : compensated 1 q = q := by simp [compensated]

/-- The positive spectral base value is computed, not assumed. All four withdrawn
majorant claims are defeated by the same admissible state at the allowed H=1. -/
theorem exceeds_every_bound (C : ℝ) :
    ∃ a : ℝ, 0 < a ∧ Admissible (velocity a) ∧
      0 < enstrophy (velocity a) ∧ 0 < palinstrophy (velocity a) ∧
      C < cauchy (velocity a) ∧
      C < compensated 1 (cauchy (velocity a)) ∧
      C < compensated 1 (data (velocity a)).schurValue ∧
      C < compensated 1 (fixedRow (velocity a) 2) := by
  obtain ⟨a, ha, _, _, hgt⟩ := inverse_exceeds (275 / 33366002112) C (by norm_num)
  have hg : C < 275 / (33366002112 * a ^ 2) := by
    convert hgt using 1
    ring
  refine ⟨a, ha, velocity_admissible a, ?_, ?_, ?_⟩
  · rw [enstrophy_velocity]; positivity
  · rw [palinstrophy_velocity]; positivity
  · simpa only [height_one, cauchy_value a ha, schur_value a ha, fixed_row_value a ha,
      and_self] using hg

/-- Unrestricted PNT-7 fails even on this fixed cutoff and threshold. -/
theorem no_cauchy_bound :
    ¬ ∃ C : ℝ, ∀ u : State, Admissible u → 0 < enstrophy u →
      0 < palinstrophy u → cauchy u ≤ C := by
  rintro ⟨C, hC⟩
  obtain ⟨a, _, hu, hZ, hP, hg, _⟩ := exceeds_every_bound C
  exact (not_lt_of_ge (hC _ hu hZ hP)) hg

/-- Unrestricted PNT-12 fails; its diagonal sum and maximum row are realized. -/
theorem no_schur_bound :
    ¬ ∃ C : ℝ, ∀ u : State, Admissible u → 0 < enstrophy u →
      0 < palinstrophy u → compensated 1 (data u).schurValue ≤ C := by
  rintro ⟨C, hC⟩
  obtain ⟨a, _, hu, hZ, hP, _, _, hg, _⟩ := exceeds_every_bound C
  exact (not_lt_of_ge (hC _ hu hZ hP)) hg

/-- Unrestricted PNT-14 fails for the actual nonzero row 2. -/
theorem no_fixed_row_bound :
    ¬ ∃ C : ℝ, ∀ u : State, Admissible u → 0 < enstrophy u →
      0 < palinstrophy u → compensated 1 (fixedRow u 2) ≤ C := by
  rintro ⟨C, hC⟩
  obtain ⟨a, _, hu, hZ, hP, _, _, _, hg⟩ := exceeds_every_bound C
  exact (not_lt_of_ge (hC _ hu hZ hP)) hg

end NavierStokesAudit.Majorants
end
