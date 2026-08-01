#pragma once

#include "local_critical_derivative_ledger.hpp"

namespace lemma {

struct LocalQuarticIdentityReport {
    SpectralReal local_advection_h1_squared = 0.0L;
    SpectralReal local_outer_state_derivative = 0.0L;
    SpectralReal local_outer_negative_square_error = 0.0L;
    SpectralReal local_enstrophy_derivative = 0.0L;
    SpectralReal local_enstrophy_identity_error = 0.0L;
    SpectralReal nonlocal_enstrophy_derivative = 0.0L;
    SpectralReal nonlocal_enstrophy_identity_error = 0.0L;
    bool finite = false;
};

class LocalQuarticIdentityLedger {
public:
    // Certifies the signed-square and enstrophy-transfer identities that
    // survive before estimating the local-local quartic SLD-1P block.
    [[nodiscard]] static LocalQuarticIdentityReport evaluate(
        const SpectralDynamics& dynamics,
        const SpectralState& state,
        const LocalCriticalDerivativeLedgerReport& derivative_ledger);
};

}  // namespace lemma
