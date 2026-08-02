#include "local_quartic_reduced_ledger.hpp"

#include <algorithm>
#include <cmath>

namespace lemma {
namespace {

SpectralReal relative_error(SpectralReal value, SpectralReal reference) {
    const SpectralReal scale = std::max(
        {std::abs(value), std::abs(reference), 1.0e-30L});
    return std::abs(value - reference) / scale;
}

}  // namespace

LocalQuarticReducedReport LocalQuarticReducedLedger::evaluate(
    const LocalQuarticCommutatorReport& commutator,
    const LocalQuarticProjectedResidualReport& projected,
    const LocalCriticalDerivativeLedgerReport& derivative_ledger) {
    LocalQuarticReducedReport result;
    result.negative_commutator_pairing =
        commutator.negative_commutator_pairing;
    result.enstrophy_normalization_remainder =
        projected.expanded_enstrophy_remainder;
    result.palinstrophy_normalization_remainder =
        projected.expanded_palinstrophy_cross;
    result.reduced_pairing = result.negative_commutator_pairing +
        result.enstrophy_normalization_remainder +
        result.palinstrophy_normalization_remainder;
    result.projected_plus_advected_pairing =
        projected.projected_pairing +
        commutator.advected_slot_derivative;
    result.raw_reconstruction_error = relative_error(
        result.reduced_pairing,
        result.projected_plus_advected_pairing);
    result.normalized_reduced_pairing =
        projected.normalized_projected_pairing +
        projected.normalized_advected_slot;
    result.normalized_advecting_slot =
        projected.normalized_advecting_slot;
    result.normalized_total = result.normalized_reduced_pairing +
        result.normalized_advecting_slot;
    result.expected_normalized_local_quartet =
        projected.expected_normalized_local_quartet;
    result.normalized_reconstruction_error = relative_error(
        result.normalized_total,
        result.expected_normalized_local_quartet);
    const SpectralReal stretching = derivative_ledger.signed_stretching;
    const SpectralReal stretching3 =
        stretching * stretching * stretching;
    const SpectralReal enstrophy = derivative_ledger.enstrophy;
    const SpectralReal palinstrophy = derivative_ledger.palinstrophy;
    const SpectralReal local_raw_bracket = result.reduced_pairing +
        derivative_ledger.stretching_local_nonlinear_roles
            .advecting_slot;
    result.polynomial_local_numerator = 4.0L * stretching3 *
        enstrophy * palinstrophy * local_raw_bracket;
    const SpectralReal enstrophy2 = enstrophy * enstrophy;
    const SpectralReal palinstrophy2 = palinstrophy * palinstrophy;
    const SpectralReal palinstrophy4 = palinstrophy2 * palinstrophy2;
    result.expected_polynomial_local_numerator =
        derivative_ledger.density_nonlinear_partition.local *
        enstrophy2 * palinstrophy4;
    result.polynomial_reconstruction_error = relative_error(
        result.polynomial_local_numerator,
        result.expected_polynomial_local_numerator);
    const SpectralReal absolute_components =
        std::abs(result.normalized_reduced_pairing) +
        std::abs(result.normalized_advecting_slot);
    if (absolute_components > 0.0L) {
        result.cancellation_fraction = 1.0L -
            std::abs(result.normalized_total) / absolute_components;
    }
    result.finite = commutator.finite && projected.finite &&
        result.raw_reconstruction_error < 1.0e-12L &&
        result.normalized_reconstruction_error < 1.0e-12L &&
        result.polynomial_reconstruction_error < 1.0e-12L;
    return result;
}

}  // namespace lemma
