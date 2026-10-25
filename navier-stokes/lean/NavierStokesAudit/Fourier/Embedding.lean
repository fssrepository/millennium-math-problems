module

public import NavierStokesAudit.Fourier.Cube
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Ring

/-! Explicit zero extension to the complete cube. The compressed six-mode input
and the full-cube energy, enstrophy, palinstrophy and pairings agree. -/
@[expose] public noncomputable section
namespace NavierStokesAudit.CubeWitness
open FourierWitness

def fullVelocity (a : ℝ) (k : Fin 27) : Vec := fun j =>
  sum6 (fun l => if supportIndex l = k then velocity a l j else 0)

theorem support_injective : Function.Injective supportIndex := by decide

theorem full_on_support (a : ℝ) (l : Fin 6) :
    fullVelocity a (supportIndex l) = velocity a l := by
  funext j
  fin_cases l <;> simp [fullVelocity, sum6, supportIndex]

theorem full_off_support (a : ℝ) (k : Fin 27) (hk : ∀ l, supportIndex l ≠ k) :
    fullVelocity a k = fun _ => 0 := by
  funext j
  simp [fullVelocity, hk, sum6]

/-- Includes reality, divergence freedom and the vanishing mean. -/
theorem full_admissible (a : ℝ) :
    (∀ k, dot (waveReal k) (fullVelocity a k) = 0) ∧
    (∀ k j, fullVelocity a k.rev j = -fullVelocity a k j) ∧
    fullVelocity a 13 = (fun _ => 0) := by
  refine ⟨?_, ?_, ?_⟩
  · intro k
    fin_cases k <;> norm_num [fullVelocity, sum6, supportIndex, dot,
      waveReal, wave, velocity, basis]
  · intro k j
    fin_cases k <;> fin_cases j <;>
      norm_num [fullVelocity, sum6, supportIndex, velocity, basis, Fin.rev]
  · funext j
    simp [fullVelocity, sum6, supportIndex]

theorem full_norms (a : ℝ) :
    sum27 (fun k => dot (fullVelocity a k) (fullVelocity a k)) = energy (velocity a) ∧
    sum27 (fun k => weight k * dot (fullVelocity a k) (fullVelocity a k)) =
      enstrophy (velocity a) ∧
    sum27 (fun k => weight k ^ 2 * dot (fullVelocity a k) (fullVelocity a k)) =
      palinstrophy (velocity a) := by
  rw [energy_velocity, enstrophy_velocity, palinstrophy_velocity]
  norm_num [sum27, fullVelocity, sum6, supportIndex, velocity, basis,
    dot, weight, lengthSq, wave]
  constructor <;> (try constructor) <;> ring

set_option maxRecDepth 4096 in
theorem full_pairings (a : ℝ) :
    sum27 (fun k => weight k * dot (fullVelocity a k)
      (advection selected (velocity a) k)) = 2 * a ^ 3 ∧
    sum27 (fun k => weight k ^ 2 * dot (fullVelocity a k)
      (advection selected (velocity a) k)) = 6 * a ^ 3 := by
  rw [advection_velocity]
  norm_num [sum27, fullVelocity, sum6, supportIndex, velocity, basis,
    advectionBasis, dot, weight, lengthSq, wave]
  constructor <;> ring

end NavierStokesAudit.CubeWitness
end
