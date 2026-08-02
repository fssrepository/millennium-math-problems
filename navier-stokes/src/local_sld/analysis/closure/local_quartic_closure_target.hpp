#pragma once

#include "local_quartic_reduced_ledger.hpp"

namespace lemma {

struct LocalQuarticClosureTargetReport {
    SpectralReal signed_two_entry_bracket = 0.0L;
    SpectralReal absolute_two_entry_bracket = 0.0L;
    SpectralReal candidate_scale = 0.0L;
    SpectralReal required_constant_ratio = 0.0L;
    SpectralReal sld_first_denominator_term = 0.0L;
    SpectralReal sld_shift_denominator_term = 0.0L;
    SpectralReal weighted_geometric_mean = 0.0L;
    SpectralReal local_absolute_polynomial_numerator = 0.0L;
    SpectralReal reconstructed_scaled_geometric_mean = 0.0L;
    SpectralReal geometric_identity_error = 0.0L;
    SpectralReal young_upper_bound = 0.0L;
    SpectralReal young_ratio = 0.0L;
    SpectralReal certified_local_constant = 0.0L;
    SpectralReal initial_frequency_power = 1.0L;
    SpectralReal initial_shift_power = 0.25L;
    SpectralReal enstrophy_power = 1.25L;
    SpectralReal palinstrophy_power = 0.75L;
    SpectralReal required_depletion_power = 0.75L;
    bool candidate_proved = false;
    bool algebra_certified = false;
    bool finite = false;
};

class LocalQuarticClosureTarget {
public:
    // Certifies the Young reduction:
    // |K+G| <= C k0 B0^(1/4) Z^(5/4)P^(3/4)
    // implies the local SLD-1P bound with A_local=3C.
    // It measures C on a finite state but does not prove a uniform C.
    [[nodiscard]] static LocalQuarticClosureTargetReport evaluate(
        const ShiftedCriticalDensityDiagnostic& diagnostic,
        const LocalCriticalDerivativeLedgerReport& derivative_ledger,
        const LocalQuarticReducedReport& reduced);
};

}  // namespace lemma
