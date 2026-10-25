module

public import NavierStokesAudit.Algebra.Normalization
public import Mathlib.Data.Fin.VecNotation
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

/-!
An exact six-mode Fourier counterexample to an amplitude-independent PNT-2 bound.

The imaginary Fourier coefficient at `wave k` is `u k`; all real parts and all
other coefficients are zero. Thus reality is antisymmetry under opposite waves.
For these imaginary-only coefficients, i (q · i u_p) (i u_q) has imaginary part
-(q · u_p) u_q. `advection` computes that convolution and the Leray projection.
All sums, norms, locality masks and the height-one core are evaluated below.
The pairings only need output modes in the velocity support; output modes
outside that support have zero velocity and contribute zero to these pairings.
-/

@[expose] public noncomputable section

namespace NavierStokesAudit.FourierWitness

abbrev Vec := Fin 3 → ℝ
abbrev State := Fin 6 → Vec

def wave : Fin 6 → Fin 3 → ℤ :=
  ![![1, 0, 0], ![0, 1, 1], ![1, 1, 1],
    ![-1, 0, 0], ![0, -1, -1], ![-1, -1, -1]]

def opposite : Fin 6 → Fin 6 := ![3, 4, 5, 0, 1, 2]

def waveReal (k : Fin 6) : Vec := fun j => (wave k j : ℝ)

def dot (u v : Vec) : ℝ := u 0 * v 0 + u 1 * v 1 + u 2 * v 2

def sum6 (f : Fin 6 → ℝ) : ℝ := f 0 + f 1 + f 2 + f 3 + f 4 + f 5

def lengthSq (k : Fin 6) : ℕ :=
  (wave k 0 ^ 2 + wave k 1 ^ 2 + wave k 2 ^ 2).toNat

def weight (k : Fin 6) : ℝ := lengthSq k

def interaction (p q k : Fin 6) : Bool :=
  (wave p 0 + wave q 0 == wave k 0) &&
  (wave p 1 + wave q 1 == wave k 1) &&
  (wave p 2 + wave q 2 == wave k 2)

def isLocal (p q k : Fin 6) : Bool :=
  decide (max (lengthSq p) (max (lengthSq q) (lengthSq k)) ≤
    4 * min (lengthSq p) (min (lengthSq q) (lengthSq k)))

/-- The maximum squared length after division by the common gcd. -/
def height (p q k : Fin 6) : ℕ :=
  max (lengthSq p) (max (lengthSq q) (lengthSq k)) /
    Nat.gcd (lengthSq p) (Nat.gcd (lengthSq q) (lengthSq k))

def localSelection (p q k : Fin 6) : Bool := interaction p q k && isLocal p q k

def coreSelection (p q k : Fin 6) : Bool :=
  localSelection p q k && decide (height p q k ≤ 1)

def rawAdvection (select : Fin 6 → Fin 6 → Fin 6 → Bool) (u : State)
    (k : Fin 6) : Vec := fun j =>
  -sum6 (fun p => sum6 (fun q =>
    if select p q k then dot (waveReal q) (u p) * u q j else 0))

def advection (select : Fin 6 → Fin 6 → Fin 6 → Bool) (u : State)
    (k : Fin 6) : Vec := fun j =>
  rawAdvection select u k j -
    dot (waveReal k) (rawAdvection select u k) / weight k * waveReal k j

def energy (u : State) : ℝ := sum6 (fun k => dot (u k) (u k))
def enstrophy (u : State) : ℝ := sum6 (fun k => weight k * dot (u k) (u k))
def palinstrophy (u : State) : ℝ := sum6 (fun k => weight k ^ 2 * dot (u k) (u k))
def stretching (u b : State) : ℝ := sum6 (fun k => weight k * dot (u k) (b k))
def crossPairing (u b : State) : ℝ := sum6 (fun k => weight k ^ 2 * dot (u k) (b k))

/-- All constraints on a real, mean-zero, divergence-free field on this support. -/
def Admissible (u : State) : Prop :=
  (∀ k, dot (waveReal k) (u k) = 0) ∧
  (∀ k j, u (opposite k) j = -u k j)

def data (u : State) : PNTData where
  energy := energy u
  Z := enstrophy u
  P := palinstrophy u
  S := stretching u (advection localSelection u)
  sc := stretching u (advection coreSelection u)
  st := stretching u (advection localSelection u) - stretching u (advection coreSelection u)
  tc := crossPairing u (advection coreSelection u)
  tt := crossPairing u (advection localSelection u) - crossPairing u (advection coreSelection u)
  -- Not used by PNT-2. No spectral realization of the Schur majorant is asserted here.
  B := 0
  D := 0
  R := 0

def basis : State :=
  ![![0, 1, 0], ![1, 0, 0], ![1, -1, 0],
    ![0, -1, 0], ![-1, 0, 0], ![-1, 1, 0]]

def velocity (a : ℝ) : State := fun k j => a * basis k j

def rawBasis : State :=
  ![![2, -1, 0], ![1, -2, 0], ![-1, -1, 0],
    ![-2, 1, 0], ![-1, 2, 0], ![1, 1, 0]]

def advectionBasis : State :=
  ![![0, -1, 0], ![1, -1, 1], ![-1/3, -1/3, 2/3],
    ![0, 1, 0], ![-1, 1, -1], ![1/3, 1/3, -2/3]]

theorem wave_geometry :
    Function.Injective wave ∧ (∀ k, 0 < lengthSq k) ∧
    (∀ k j, wave (opposite k) j = -wave k j) := by
  decide

theorem interaction_geometry :
    ∀ p q k, interaction p q k = true → isLocal p q k = true ∧ height p q k = 3 := by
  decide

theorem core_empty : ∀ p q k, coreSelection p q k = false := by
  decide

theorem velocity_admissible (a : ℝ) : Admissible (velocity a) := by
  constructor
  · intro k
    fin_cases k <;> simp [dot, waveReal, wave, velocity, basis]
  · intro k j
    fin_cases k <;> fin_cases j <;> simp [velocity, opposite, basis]

set_option maxRecDepth 4096 in
set_option maxHeartbeats 2000000 in
theorem raw_velocity (a : ℝ) :
    rawAdvection localSelection (velocity a) = fun k j => a ^ 2 * rawBasis k j := by
  funext k j
  fin_cases k <;> fin_cases j <;>
    norm_num [rawAdvection, sum6, localSelection, interaction, isLocal, lengthSq,
      wave, dot, waveReal, velocity, basis, rawBasis] <;> ring

theorem advection_velocity (a : ℝ) :
    advection localSelection (velocity a) = fun k j => a ^ 2 * advectionBasis k j := by
  funext k j
  simp only [advection, raw_velocity]
  fin_cases k <;> fin_cases j <;>
    norm_num [dot, waveReal, weight, lengthSq, wave, rawBasis, advectionBasis] <;> ring

theorem core_advection (u : State) : advection coreSelection u = fun _ _ => 0 := by
  funext k j
  simp [advection, rawAdvection, core_empty, sum6, dot]

theorem energy_velocity (a : ℝ) : energy (velocity a) = 8 * a ^ 2 := by
  norm_num [energy, sum6, dot, velocity, basis]
  ring

theorem enstrophy_velocity (a : ℝ) : enstrophy (velocity a) = 18 * a ^ 2 := by
  norm_num [enstrophy, sum6, dot, velocity, basis, weight, lengthSq, wave]
  ring

theorem palinstrophy_velocity (a : ℝ) : palinstrophy (velocity a) = 46 * a ^ 2 := by
  norm_num [palinstrophy, sum6, dot, velocity, basis, weight, lengthSq, wave]
  ring

theorem stretching_velocity (a : ℝ) :
    stretching (velocity a) (advection localSelection (velocity a)) = 2 * a ^ 3 := by
  rw [advection_velocity]
  norm_num [stretching, sum6, dot, velocity, basis, advectionBasis, weight, lengthSq, wave]
  ring

theorem cross_velocity (a : ℝ) :
    crossPairing (velocity a) (advection localSelection (velocity a)) = 6 * a ^ 3 := by
  rw [advection_velocity]
  norm_num [crossPairing, sum6, dot, velocity, basis, advectionBasis, weight, lengthSq, wave]
  ring

def baseData : PNTData := ⟨8, 18, 46, 2, 0, 2, 0, 6, 0, 0, 0⟩

theorem data_velocity (a : ℝ) : data (velocity a) = baseData.scale a := by
  unfold data
  rw [energy_velocity, enstrophy_velocity, palinstrophy_velocity,
    stretching_velocity, cross_velocity, core_advection]
  simp [stretching, crossPairing, sum6, dot, baseData, PNTData.scale, mul_comm]

theorem base_value : baseData.openValue = 1 / 876024 := by
  norm_num [baseData, PNTData.openValue, PNTData.product]

/-- Exact value of the repository's PNT-2 formula on an admissible amplitude ray. -/
theorem value (a : ℝ) (ha : 0 < a) :
    (data (velocity a)).openValue = 1 / (876024 * a) := by
  rw [data_velocity, open_scale baseData a ha (by norm_num [baseData])
    (by norm_num [baseData]), base_value]
  ring

/-- Every proposed real bound is exceeded by a member of this explicit family. -/
theorem exceeds_every_bound (C : ℝ) :
    ∃ a : ℝ, 0 < a ∧ Admissible (velocity a) ∧
      0 < enstrophy (velocity a) ∧ 0 < palinstrophy (velocity a) ∧
      C < (data (velocity a)).openValue := by
  obtain ⟨a, ha, _, hgt, _⟩ := inverse_exceeds (1 / 876024) C (by norm_num)
  refine ⟨a, ha, velocity_admissible a, ?_, ?_, ?_⟩
  · rw [enstrophy_velocity]; positivity
  · rw [palinstrophy_velocity]; positivity
  · rw [value a ha]
    convert hgt using 1
    ring

/-- The unrestricted PNT-2 bound is false already on this six-mode support. -/
theorem no_uniform_bound :
    ¬ ∃ C : ℝ, ∀ u : State, Admissible u → 0 < enstrophy u →
      0 < palinstrophy u → (data u).openValue ≤ C := by
  rintro ⟨C, hC⟩
  obtain ⟨a, _, hu, hZ, hP, hgt⟩ := exceeds_every_bound C
  exact (not_lt_of_ge (hC (velocity a) hu hZ hP)) hgt

end NavierStokesAudit.FourierWitness

end
