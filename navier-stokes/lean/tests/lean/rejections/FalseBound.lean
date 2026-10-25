module

import NavierStokesAudit
import Mathlib.Tactic.NormNum

open NavierStokesAudit.FourierWitness

-- Expected failure: halving the amplitude doubles the open quotient.
example : (data (velocity (1/2))).openValue ≤ (data (velocity 1)).openValue := by
  rw [value (1/2) (by norm_num), value 1 (by norm_num)]
  norm_num
