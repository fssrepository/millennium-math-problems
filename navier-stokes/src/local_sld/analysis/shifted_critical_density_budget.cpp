#include "shifted_critical_density_budget.hpp"

#include <algorithm>
#include <cmath>

namespace lemma {

ShiftedCriticalDensityBudget
ShiftedCriticalDensityBudgetAnalyzer::evaluate(
    const ShiftedCriticalDensityDiagnostic& diagnostic,
    const LocalCriticalDerivativeLedgerReport& ledger) {
    ShiftedCriticalDensityBudget result;
    result.shifted_density =
        diagnostic.local_critical_density + diagnostic.initial_ep_shift;
    if (!(result.shifted_density > 0.0L) ||
        !(diagnostic.normalization > 0.0L)) {
        return result;
    }
    result.density_fraction =
        diagnostic.local_critical_density / result.shifted_density;
    result.log_rate_from_stretching =
        ledger.density_from_stretching.total / result.shifted_density;
    result.log_rate_from_enstrophy =
        ledger.density_from_enstrophy.total / result.shifted_density;
    result.log_rate_from_palinstrophy =
        ledger.density_from_palinstrophy.total / result.shifted_density;
    result.normalized_from_stretching =
        result.log_rate_from_stretching / diagnostic.normalization;
    result.normalized_from_enstrophy =
        result.log_rate_from_enstrophy / diagnostic.normalization;
    result.normalized_from_palinstrophy =
        result.log_rate_from_palinstrophy / diagnostic.normalization;
    result.normalized_nonlinear =
        ledger.reconstructed_density_derivative.nonlinear /
        (result.shifted_density * diagnostic.normalization);
    result.normalized_local_nonlinear =
        ledger.density_nonlinear_partition.local /
        (result.shifted_density * diagnostic.normalization);
    result.normalized_nonlocal_nonlinear =
        ledger.density_nonlinear_partition.nonlocal /
        (result.shifted_density * diagnostic.normalization);
    if (ledger.enstrophy > 0.0L && ledger.palinstrophy > 0.0L) {
        const SpectralReal stretching2 =
            ledger.signed_stretching * ledger.signed_stretching;
        const SpectralReal stretching3 =
            stretching2 * ledger.signed_stretching;
        const SpectralReal stretching4 = stretching2 * stretching2;
        const SpectralReal palinstrophy3 =
            ledger.palinstrophy * ledger.palinstrophy *
            ledger.palinstrophy;
        const SpectralReal scale =
            result.shifted_density * diagnostic.normalization;
        const SpectralReal stretching_coefficient =
            4.0L * stretching3 /
            (ledger.enstrophy * palinstrophy3 * scale);
        const SpectralReal enstrophy_coefficient =
            -stretching4 /
            (ledger.enstrophy * ledger.enstrophy *
             palinstrophy3 * scale);
        const SpectralReal palinstrophy_coefficient =
            -3.0L * stretching4 /
            (ledger.enstrophy * palinstrophy3 *
             ledger.palinstrophy * scale);
        result.normalized_local_outer_state = stretching_coefficient *
            ledger.stretching_local_nonlinear_roles.outer_state;
        result.normalized_local_advecting_slot = stretching_coefficient *
            ledger.stretching_local_nonlinear_roles.advecting_slot;
        result.normalized_local_advected_slot = stretching_coefficient *
            ledger.stretching_local_nonlinear_roles.advected_slot;
        result.normalized_local_enstrophy = enstrophy_coefficient *
            ledger.enstrophy_nonlinear_partition.local;
        result.normalized_local_palinstrophy = palinstrophy_coefficient *
            ledger.palinstrophy_nonlinear_partition.local;
    }
    result.reconstructed_local_nonlinear =
        result.normalized_local_outer_state +
        result.normalized_local_advecting_slot +
        result.normalized_local_advected_slot +
        result.normalized_local_enstrophy +
        result.normalized_local_palinstrophy;
    const SpectralReal local_scale = std::max(
        {std::abs(result.reconstructed_local_nonlinear),
         std::abs(result.normalized_local_nonlinear), 1.0e-30L});
    result.local_nonlinear_reconstruction_error = std::abs(
        result.reconstructed_local_nonlinear -
        result.normalized_local_nonlinear) / local_scale;
    result.normalized_viscous =
        ledger.reconstructed_density_derivative.viscous /
        (result.shifted_density * diagnostic.normalization);
    result.reconstructed_normalized_rate =
        result.normalized_from_stretching +
        result.normalized_from_enstrophy +
        result.normalized_from_palinstrophy;
    const SpectralReal reference =
        diagnostic.normalized_shifted_log_derivative;
    const SpectralReal scale = std::max(
        {std::abs(result.reconstructed_normalized_rate),
         std::abs(reference), 1.0e-30L});
    result.relative_reconstruction_error = std::abs(
        result.reconstructed_normalized_rate - reference) / scale;
    const SpectralReal stretching2 =
        ledger.signed_stretching * ledger.signed_stretching;
    const SpectralReal stretching3 =
        stretching2 * ledger.signed_stretching;
    const SpectralReal stretching4 = stretching2 * stretching2;
    result.polynomial_numerator =
        4.0L * stretching3 * ledger.stretching_derivative.total *
            ledger.enstrophy * ledger.palinstrophy -
        stretching4 * ledger.enstrophy_derivative.total *
            ledger.palinstrophy -
        3.0L * stretching4 * ledger.enstrophy *
            ledger.palinstrophy_derivative.total;
    result.polynomial_denominator = diagnostic.initial_frequency *
        (stretching4 * ledger.enstrophy * ledger.enstrophy *
             ledger.palinstrophy +
         diagnostic.initial_ep_shift * ledger.enstrophy *
             ledger.enstrophy * ledger.enstrophy *
             ledger.palinstrophy * ledger.palinstrophy *
             ledger.palinstrophy * ledger.palinstrophy);
    if (result.polynomial_denominator > 0.0L) {
        result.polynomial_required_coefficient =
            result.polynomial_numerator /
            result.polynomial_denominator;
    }
    const SpectralReal polynomial_scale = std::max(
        {std::abs(result.polynomial_required_coefficient),
         std::abs(reference), 1.0e-30L});
    result.polynomial_equivalence_error = std::abs(
        result.polynomial_required_coefficient - reference) /
        polynomial_scale;
    result.finite =
        ledger.finite &&
        std::isfinite(result.reconstructed_normalized_rate) &&
        result.relative_reconstruction_error < 1.0e-12L &&
        result.local_nonlinear_reconstruction_error < 1.0e-12L &&
        result.polynomial_equivalence_error < 1.0e-12L;
    return result;
}

}  // namespace lemma
