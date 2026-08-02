#pragma once

#include "local_quartic_commutator.hpp"
#include "local_quartic_projected_residual.hpp"

namespace lemma {

struct LocalQuarticReducedReport {
    SpectralReal negative_commutator_pairing = 0.0L;
    SpectralReal enstrophy_normalization_remainder = 0.0L;
    SpectralReal palinstrophy_normalization_remainder = 0.0L;
    SpectralReal reduced_pairing = 0.0L;
    SpectralReal projected_plus_advected_pairing = 0.0L;
    SpectralReal raw_reconstruction_error = 0.0L;
    SpectralReal normalized_reduced_pairing = 0.0L;
    SpectralReal normalized_advecting_slot = 0.0L;
    SpectralReal normalized_total = 0.0L;
    SpectralReal expected_normalized_local_quartet = 0.0L;
    SpectralReal normalized_reconstruction_error = 0.0L;
    SpectralReal polynomial_local_numerator = 0.0L;
    SpectralReal expected_polynomial_local_numerator = 0.0L;
    SpectralReal polynomial_reconstruction_error = 0.0L;
    SpectralReal cancellation_fraction = 0.0L;
    bool finite = false;
};

class LocalQuarticReducedLedger {
public:
    // Combines the commutator and normalization identities so that the
    // complete local quartet has only two entries: the reduced pairing and
    // the advecting Frechet slot.
    [[nodiscard]] static LocalQuarticReducedReport evaluate(
        const LocalQuarticCommutatorReport& commutator,
        const LocalQuarticProjectedResidualReport& projected,
        const LocalCriticalDerivativeLedgerReport& derivative_ledger);
};

}  // namespace lemma
