module

import NavierStokesAudit
import ProofAudit

namespace NavierStokesAudit

axiom injected : False
theorem injected_result : False := injected

end NavierStokesAudit

-- Expected failure even though the preceding theorem typechecks.
audit_pnt_proofs
