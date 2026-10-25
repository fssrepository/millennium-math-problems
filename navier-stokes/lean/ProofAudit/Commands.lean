module

public meta import Lean.Util.CollectAxioms
public meta import Lean.Elab.Command

public meta section

open Lean Elab Command

/-- Audit every declaration in our namespace, including transitive dependencies.
Use Lean's exported axiom metadata: the module system can erase imported proof
bodies and present their declarations as axioms in the elaborator environment.
Those placeholders are not new trusted axioms. -/
elab "audit_pnt_proofs" : command => do
  let env ← getEnv
  let allowed : Array Name := #[``propext, ``Classical.choice, ``Quot.sound]
  let mut count : Nat := 0
  for (name, _) in env.constants.toList do
    if (`NavierStokesAudit).isPrefixOf (privateToUserName name) then
      count := count + 1
      for axiomName in ← collectAxioms name do
        unless allowed.contains axiomName do
          throwError "PNT audit: unexpected axiom {axiomName} in {name}"
  if count == 0 then
    throwError "PNT audit: no project declarations were loaded"
  for required in #[`NavierStokesAudit.open_unbounded,
      `NavierStokesAudit.schur_unbounded,
      `NavierStokesAudit.FourierWitness.value,
      `NavierStokesAudit.FourierWitness.no_uniform_bound,
      `NavierStokesAudit.CubeWitness.full_admissible,
      `NavierStokesAudit.CubeWitness.full_norms,
      `NavierStokesAudit.CubeWitness.full_pairings,
      `NavierStokesAudit.Cauchy.factorization,
      `NavierStokesAudit.Cauchy.majorant_square,
      `NavierStokesAudit.Schur.vector_bound,
      `NavierStokesAudit.Majorants.exceeds_every_bound,
      `NavierStokesAudit.Majorants.no_cauchy_bound,
      `NavierStokesAudit.Majorants.no_schur_bound,
      `NavierStokesAudit.Majorants.no_fixed_row_bound] do
    unless (env.find? required).isSome do
      throwError "PNT audit: missing required theorem {required}"
  logInfo m!"PNT audit PASS: {count} declarations; only propext, Classical.choice, Quot.sound allowed"

end
