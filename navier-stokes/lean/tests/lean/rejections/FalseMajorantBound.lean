module

import NavierStokesAudit
import Mathlib.Tactic.NormNum

open NavierStokesAudit.FourierWitness
open NavierStokesAudit.Majorants

-- Expected failure: halving the amplitude multiplies the squared majorant by four.
example : cauchy (velocity (1/2)) ≤ cauchy (velocity 1) := by
  rw [cauchy_value (1/2) (by norm_num), cauchy_value 1 (by norm_num)]
  norm_num
