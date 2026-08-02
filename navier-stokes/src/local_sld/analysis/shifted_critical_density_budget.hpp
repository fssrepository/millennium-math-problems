#pragma once

#include "local_critical_derivative_ledger.hpp"
#include "shifted_critical_density.hpp"

namespace lemma {

struct ShiftedCriticalDensityBudget {
    SpectralReal shifted_density = 0.0L;
    SpectralReal density_fraction = 0.0L;
    SpectralReal log_rate_from_stretching = 0.0L;
    SpectralReal log_rate_from_enstrophy = 0.0L;
    SpectralReal log_rate_from_palinstrophy = 0.0L;
    SpectralReal normalized_from_stretching = 0.0L;
    SpectralReal normalized_from_enstrophy = 0.0L;
    SpectralReal normalized_from_palinstrophy = 0.0L;
    SpectralReal normalized_nonlinear = 0.0L;
    SpectralReal normalized_local_nonlinear = 0.0L;
    SpectralReal normalized_nonlocal_nonlinear = 0.0L;
    SpectralReal normalized_local_outer_state = 0.0L;
    SpectralReal normalized_local_advecting_slot = 0.0L;
    SpectralReal normalized_local_advected_slot = 0.0L;
    SpectralReal normalized_local_enstrophy = 0.0L;
    SpectralReal normalized_local_palinstrophy = 0.0L;
    SpectralReal reconstructed_local_nonlinear = 0.0L;
    SpectralReal local_nonlinear_reconstruction_error = 0.0L;
    SpectralReal normalized_viscous = 0.0L;
    SpectralReal reconstructed_normalized_rate = 0.0L;
    SpectralReal relative_reconstruction_error = 0.0L;
    SpectralReal polynomial_numerator = 0.0L;
    SpectralReal polynomial_denominator = 0.0L;
    SpectralReal polynomial_required_coefficient = 0.0L;
    SpectralReal polynomial_equivalence_error = 0.0L;
    bool finite = false;
};

class ShiftedCriticalDensityBudgetAnalyzer {
public:
    // Computes the exact dimensionless SLD-1 budget. Using chain-rule
    // density contributions instead of S'/S keeps the formula regular at
    // the zero-stretching set where C_local=0.
    [[nodiscard]] static ShiftedCriticalDensityBudget evaluate(
        const ShiftedCriticalDensityDiagnostic& diagnostic,
        const LocalCriticalDerivativeLedgerReport& ledger);
};

}  // namespace lemma
