#include "local_quartic_projected_residual.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace lemma {
namespace {

SpectralReal relative_error(SpectralReal value, SpectralReal reference) {
    const SpectralReal scale = std::max(
        {std::abs(value), std::abs(reference), 1.0e-30L});
    return std::abs(value - reference) / scale;
}

}  // namespace

LocalQuarticProjectedResidualReport
LocalQuarticProjectedResidual::evaluate(
    const SpectralDynamics& dynamics,
    const SpectralState& state,
    const ShiftedCriticalDensityDiagnostic& diagnostic,
    const LocalCriticalDerivativeLedgerReport& derivative_ledger) {
    LocalQuarticProjectedResidualReport result;
    const SpectralReal stretching = derivative_ledger.signed_stretching;
    const SpectralReal enstrophy = derivative_ledger.enstrophy;
    const SpectralReal palinstrophy = derivative_ledger.palinstrophy;
    if (!(enstrophy > 0.0L && palinstrophy > 0.0L)) {
        return result;
    }
    const SpectralIncrement local_advection =
        dynamics.advection_direct_partition(
            state, TriadPartition::local);
    const SpectralReal z_coefficient = stretching / (2.0L * enstrophy);
    const SpectralReal p_coefficient =
        3.0L * stretching / (2.0L * palinstrophy);
    result.completion_coefficient =
        3.0L * stretching / (4.0L * palinstrophy);

    SpectralReal negative_square = 0.0L;
    SpectralReal palinstrophy_cross_base = 0.0L;
    SpectralReal completed_square = 0.0L;
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const SpectralReal wave2 = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        const SpectralReal wave4 = wave2 * wave2;
        ComplexVector residual{};
        ComplexVector completed_residual{};
        for (std::size_t component = 0; component < 3; ++component) {
            residual[component] =
                local_advection[mode][component] -
                z_coefficient * state.velocity[mode][component] -
                p_coefficient * wave2 *
                    state.velocity[mode][component];
            completed_residual[component] =
                local_advection[mode][component] -
                result.completion_coefficient * wave2 *
                    state.velocity[mode][component];
        }
        result.projected_pairing -= wave2 * std::real(dot_hermitian(
            local_advection[mode], residual));
        negative_square -= wave2 * std::real(dot_hermitian(
            local_advection[mode], local_advection[mode]));
        palinstrophy_cross_base += wave4 * std::real(dot_hermitian(
            local_advection[mode], state.velocity[mode]));
        completed_square -= wave2 * std::real(dot_hermitian(
            completed_residual, completed_residual));
    }
    result.expanded_negative_square = negative_square;
    result.expanded_enstrophy_remainder =
        stretching * stretching / (2.0L * enstrophy);
    result.expanded_palinstrophy_cross =
        p_coefficient * palinstrophy_cross_base;
    result.expanded_total = result.expanded_negative_square +
        result.expanded_enstrophy_remainder +
        result.expanded_palinstrophy_cross;
    result.completed_negative_square = completed_square;
    result.completed_enstrophy_remainder =
        result.expanded_enstrophy_remainder;
    result.completed_hyperpalinstrophy_remainder =
        result.completion_coefficient * result.completion_coefficient *
        derivative_ledger.hyperpalinstrophy;
    result.completed_total = result.completed_negative_square +
        result.completed_enstrophy_remainder +
        result.completed_hyperpalinstrophy_remainder;
    result.expansion_relative_error = relative_error(
        result.projected_pairing, result.expanded_total);
    result.completion_relative_error = relative_error(
        result.projected_pairing, result.completed_total);

    const SpectralReal shifted_density =
        diagnostic.local_critical_density + diagnostic.initial_ep_shift;
    const SpectralReal normalization_scale =
        shifted_density * diagnostic.normalization;
    if (normalization_scale > 0.0L) {
        const SpectralReal stretching2 = stretching * stretching;
        const SpectralReal stretching3 = stretching2 * stretching;
        const SpectralReal palinstrophy3 =
            palinstrophy * palinstrophy * palinstrophy;
        const SpectralReal normalized_coefficient =
            4.0L * stretching3 /
            (enstrophy * palinstrophy3 * normalization_scale);
        result.normalized_projected_pairing =
            normalized_coefficient * result.projected_pairing;
        result.normalized_advecting_slot = normalized_coefficient *
            derivative_ledger.stretching_local_nonlinear_roles
                .advecting_slot;
        result.normalized_advected_slot = normalized_coefficient *
            derivative_ledger.stretching_local_nonlinear_roles
                .advected_slot;
        result.normalized_total =
            result.normalized_projected_pairing +
            result.normalized_advecting_slot +
            result.normalized_advected_slot;
        result.expected_normalized_local_quartet =
            derivative_ledger.density_nonlinear_partition.local /
            normalization_scale;
        result.normalized_reconstruction_error = relative_error(
            result.normalized_total,
            result.expected_normalized_local_quartet);
    }
    result.finite =
        std::isfinite(result.projected_pairing) &&
        result.expansion_relative_error < 1.0e-12L &&
        result.completion_relative_error < 1.0e-12L &&
        result.normalized_reconstruction_error < 1.0e-12L;
    return result;
}

}  // namespace lemma
