#pragma once

#include "local_critical_derivative_ledger.hpp"
#include "shifted_critical_density.hpp"

namespace lemma {

struct LocalQuarticProjectedResidualReport {
    SpectralReal projected_pairing = 0.0L;
    SpectralReal expanded_negative_square = 0.0L;
    SpectralReal expanded_enstrophy_remainder = 0.0L;
    SpectralReal expanded_palinstrophy_cross = 0.0L;
    SpectralReal expanded_total = 0.0L;
    SpectralReal completed_negative_square = 0.0L;
    SpectralReal completed_enstrophy_remainder = 0.0L;
    SpectralReal completed_hyperpalinstrophy_remainder = 0.0L;
    SpectralReal completed_total = 0.0L;
    SpectralReal completion_coefficient = 0.0L;
    SpectralReal expansion_relative_error = 0.0L;
    SpectralReal completion_relative_error = 0.0L;
    SpectralReal normalized_projected_pairing = 0.0L;
    SpectralReal normalized_advecting_slot = 0.0L;
    SpectralReal normalized_advected_slot = 0.0L;
    SpectralReal normalized_total = 0.0L;
    SpectralReal expected_normalized_local_quartet = 0.0L;
    SpectralReal normalized_reconstruction_error = 0.0L;
    bool finite = false;
};

class LocalQuarticProjectedResidual {
public:
    // For h=-B_local, combines the outer, Z', and P' chain-rule terms:
    // <A h, B-Su/(2Z)-3S Au/(2P)>.
    // It also certifies the corresponding completed-square identity.
    [[nodiscard]] static LocalQuarticProjectedResidualReport evaluate(
        const SpectralDynamics& dynamics,
        const SpectralState& state,
        const ShiftedCriticalDensityDiagnostic& diagnostic,
        const LocalCriticalDerivativeLedgerReport& derivative_ledger);
};

}  // namespace lemma
