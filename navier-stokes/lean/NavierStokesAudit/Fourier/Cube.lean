module

public import NavierStokesAudit.Fourier.SixMode
public import Mathlib.Data.Fintype.Fin
import Mathlib.Data.Fintype.Pi
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.IntervalCases
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

/-! Complete cutoff-one output cube. Inputs use the six-mode support; every
other input coefficient is zero. The zero output is explicitly masked. The
convolution therefore includes every nonzero input pair and every retained
output, including the two outputs outside the velocity support. -/
@[expose] public noncomputable section
namespace NavierStokesAudit.CubeWitness
open FourierWitness

abbrev Output := Fin 27 → Vec

def wave : Fin 27 → Fin 3 → ℤ :=
    ![![-1, -1, -1],
    ![-1, -1, 0],
    ![-1, -1, 1],
    ![-1, 0, -1],
    ![-1, 0, 0],
    ![-1, 0, 1],
    ![-1, 1, -1],
    ![-1, 1, 0],
    ![-1, 1, 1],
    ![0, -1, -1],
    ![0, -1, 0],
    ![0, -1, 1],
    ![0, 0, -1],
    ![0, 0, 0],
    ![0, 0, 1],
    ![0, 1, -1],
    ![0, 1, 0],
    ![0, 1, 1],
    ![1, -1, -1],
    ![1, -1, 0],
    ![1, -1, 1],
    ![1, 0, -1],
    ![1, 0, 0],
    ![1, 0, 1],
    ![1, 1, -1],
    ![1, 1, 0],
    ![1, 1, 1]]

def lengthSq (k : Fin 27) : ℕ :=
  (wave k 0 ^ 2 + wave k 1 ^ 2 + wave k 2 ^ 2).toNat

def weight (k : Fin 27) : ℝ := lengthSq k

def waveReal (k : Fin 27) : Vec := fun j => (wave k j : ℝ)

def interaction (p q : Fin 6) (k : Fin 27) : Bool :=
  decide (0 < lengthSq k) &&
  (FourierWitness.wave p 0 + FourierWitness.wave q 0 == wave k 0) &&
  (FourierWitness.wave p 1 + FourierWitness.wave q 1 == wave k 1) &&
  (FourierWitness.wave p 2 + FourierWitness.wave q 2 == wave k 2)

def isLocal (p q : Fin 6) (k : Fin 27) : Bool :=
  decide (max (FourierWitness.lengthSq p) (max (FourierWitness.lengthSq q) (lengthSq k)) ≤
    4 * min (FourierWitness.lengthSq p) (min (FourierWitness.lengthSq q) (lengthSq k)))

def height (p q : Fin 6) (k : Fin 27) : ℕ :=
  max (FourierWitness.lengthSq p) (max (FourierWitness.lengthSq q) (lengthSq k)) /
    Nat.gcd (FourierWitness.lengthSq p) (Nat.gcd (FourierWitness.lengthSq q) (lengthSq k))

def selected (p q : Fin 6) (k : Fin 27) : Bool := interaction p q k && isLocal p q k

/-- Dyadic height shell j has 2^(j-1) < height ≤ 2^j for j > 0. -/
def inShell (j h : ℕ) : Bool :=
  if j = 0 then decide (h ≤ 1) else decide (2 ^ (j-1) < h ∧ h ≤ 2 ^ j)

def shellSelection (j : ℕ) (p q : Fin 6) (k : Fin 27) : Bool :=
  selected p q k && decide (1 < height p q k) && inShell j (height p q k)

def rawAdvection (select : Fin 6 → Fin 6 → Fin 27 → Bool) (u : State)
    (k : Fin 27) : Vec := fun j =>
  -sum6 (fun p => sum6 (fun q =>
    if select p q k then dot (FourierWitness.waveReal q) (u p) * u q j else 0))

def advection (select : Fin 6 → Fin 6 → Fin 27 → Bool) (u : State)
    (k : Fin 27) : Vec := fun j =>
  rawAdvection select u k j -
    dot (waveReal k) (rawAdvection select u k) / weight k * waveReal k j

def sum27 (f : Fin 27 → ℝ) : ℝ := f 0 + f 1 + f 2 + f 3 + f 4 + f 5 + f 6 + f 7 + f 8 + f 9 + f 10 + f 11 + f 12 + f 13 + f 14 + f 15 + f 16 + f 17 + f 18 + f 19 + f 20 + f 21 + f 22 + f 23 + f 24 + f 25 + f 26

def h1 (b : Output) : ℝ := sum27 (fun k => weight k * dot (b k) (b k))
def gram (b c : Output) : ℝ := sum27 (fun k => weight k ^ 2 * dot (b k) (c k))
def h2 (b : Output) : ℝ := gram b b

def rawBasis : Output := ![![1, 1, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![-2, 1, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![1, -1, 0],
    ![-1, 2, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![1, -2, 0],
    ![-1, 1, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![2, -1, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![-1, -1, 0]]

def advectionBasis : Output := ![![1/3, 1/3, -2/3],
    ![0, 0, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![0, 1, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![1/3, -1/3, 2/3],
    ![-1, 1, -1],
    ![0, 0, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![1, -1, 1],
    ![-1/3, 1/3, -2/3],
    ![0, 0, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![0, -1, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![0, 0, 0],
    ![-1/3, -1/3, 2/3]]

/-- No retained nonzero cube wave is omitted or duplicated. -/
theorem cube_complete : Function.Injective wave ∧
    ∀ x y z : Fin 3, ∃ k, wave k = ![(x.val : ℤ)-1, (y.val : ℤ)-1, (z.val : ℤ)-1] := by
  decide

set_option maxRecDepth 8192 in
set_option maxHeartbeats 4000000 in
theorem interaction_geometry : ∀ p q k, interaction p q k = true →
    isLocal p q k = true ∧ height p q k = 3 := by
  decide

theorem shell_two : shellSelection 2 = selected := by
  funext p q k
  unfold shellSelection selected
  cases hi : interaction p q k <;> simp
  have hg := interaction_geometry p q k hi
  simp [hg, inShell]

theorem height_three_shell (j : ℕ) : inShell j 3 = true ↔ j = 2 := by
  by_cases hj : j = 0
  · subst j; norm_num [inShell]
  simp only [inShell, hj, ↓reduceIte, decide_eq_true_eq]
  constructor
  · intro h
    have hp : j-1 < 2 ^ (j-1) := Nat.lt_two_pow_self
    have hj3 : j ≤ 3 := by omega
    interval_cases j <;> norm_num at *
  · rintro rfl
    norm_num

theorem other_shell_empty (j : ℕ) (hj : j ≠ 2) :
    shellSelection j = fun _ _ _ => false := by
  funext p q k
  unfold shellSelection selected
  cases hi : interaction p q k <;> simp
  have hg := interaction_geometry p q k hi
  have hh : inShell j 3 = false := by
    cases h : inShell j 3
    · rfl
    · exact False.elim (hj ((height_three_shell j).mp h))
  simp [hg, hh]

theorem other_shell_advection (j : ℕ) (hj : j ≠ 2) (u : State) :
    advection (shellSelection j) u = fun _ _ => 0 := by
  rw [other_shell_empty j hj]
  funext k l
  simp [advection, rawAdvection, sum6, dot]

set_option maxRecDepth 8192 in
set_option maxHeartbeats 4000000 in
theorem raw_velocity (a : ℝ) :
    rawAdvection selected (velocity a) = fun k j => a ^ 2 * rawBasis k j := by
  funext k j
  fin_cases k <;> fin_cases j <;>
    norm_num [rawAdvection, sum6, selected, interaction, isLocal, lengthSq,
      wave, FourierWitness.lengthSq, FourierWitness.wave, dot,
      FourierWitness.waveReal, velocity, basis, rawBasis] <;> ring

set_option maxRecDepth 8192 in
set_option maxHeartbeats 2000000 in
theorem advection_velocity (a : ℝ) :
    advection selected (velocity a) = fun k j => a ^ 2 * advectionBasis k j := by
  funext k j
  simp only [advection, raw_velocity]
  fin_cases k <;> fin_cases j <;>
    norm_num [dot, waveReal, weight, lengthSq, wave, rawBasis, advectionBasis] <;> ring

theorem h1_velocity (a : ℝ) : h1 (advection selected (velocity a)) = 22 * a ^ 4 := by
  rw [advection_velocity]
  norm_num [h1, sum27, dot, weight, lengthSq, wave, advectionBasis]
  ring

theorem h2_velocity (a : ℝ) : h2 (advection selected (velocity a)) = 50 * a ^ 4 := by
  rw [advection_velocity]
  norm_num [h2, gram, sum27, dot, weight, lengthSq, wave, advectionBasis]
  ring

def supportIndex : Fin 6 → Fin 27 := ![22, 17, 26, 4, 9, 0]

theorem support_waves : ∀ k, wave (supportIndex k) = FourierWitness.wave k := by
  decide

theorem support_advection (a : ℝ) (k : Fin 6) :
    advection selected (velocity a) (supportIndex k) =
      FourierWitness.advection localSelection (velocity a) k := by
  rw [advection_velocity, FourierWitness.advection_velocity]
  funext j
  fin_cases k <;> fin_cases j <;> norm_num [supportIndex, advectionBasis,
    FourierWitness.advectionBasis]

end NavierStokesAudit.CubeWitness
end
