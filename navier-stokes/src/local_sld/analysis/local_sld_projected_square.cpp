#include "local_sld_projected_square.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace lemma {
namespace {

SpectralReal relative_error(SpectralReal left, SpectralReal right) {
    const SpectralReal scale = std::max(
        {std::abs(left), std::abs(right), 1.0e-30L});
    return std::abs(left - right) / scale;
}

}  // namespace

LocalSldProjectedSquareReport LocalSldProjectedSquare::evaluate(
    const SpectralDynamics& dynamics,
    const SpectralState& state,
    TriadSelection selection) {
    LocalSldProjectedSquareReport result;
    const SpectralIncrement advection =
        dynamics.advection_direct_partition(state, selection);

    SpectralReal advection_h1_norm2 = 0.0L;
    SpectralReal completed_h1_norm2 = 0.0L;
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const SpectralReal wave2 = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        const SpectralReal wave4 = wave2 * wave2;
        const SpectralReal wave6 = wave4 * wave2;
        const SpectralReal velocity2 = std::real(dot_hermitian(
            state.velocity[mode], state.velocity[mode]));
        const SpectralReal advection_velocity = std::real(
            dot_hermitian(advection[mode], state.velocity[mode]));
        const SpectralReal advection2 = std::real(dot_hermitian(
            advection[mode], advection[mode]));
        result.enstrophy += wave2 * velocity2;
        result.palinstrophy += wave4 * velocity2;
        result.hyperpalinstrophy += wave6 * velocity2;
        result.stretching += wave2 * advection_velocity;
        result.palinstrophy_cross += wave4 * advection_velocity;
        advection_h1_norm2 += wave2 * advection2;
    }
    if (!(result.enstrophy > 0.0L && result.palinstrophy > 0.0L)) {
        return result;
    }

    result.completion_coefficient =
        3.0L * result.stretching / (4.0L * result.palinstrophy);
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const SpectralReal wave2 = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        ComplexVector residual{};
        for (std::size_t component = 0; component < 3; ++component) {
            residual[component] = advection[mode][component] -
                result.completion_coefficient * wave2 *
                    state.velocity[mode][component];
        }
        completed_h1_norm2 += wave2 * std::real(dot_hermitian(
            residual, residual));
    }

    result.expanded_negative_square = -advection_h1_norm2;
    result.expanded_enstrophy_remainder =
        result.stretching * result.stretching /
        (2.0L * result.enstrophy);
    result.expanded_palinstrophy_remainder =
        3.0L * result.stretching * result.palinstrophy_cross /
        (2.0L * result.palinstrophy);
    result.expanded_total = result.expanded_negative_square +
        result.expanded_enstrophy_remainder +
        result.expanded_palinstrophy_remainder;

    result.completed_negative_square = -completed_h1_norm2;
    result.completed_enstrophy_remainder =
        result.expanded_enstrophy_remainder;
    result.completed_hyperpalinstrophy_remainder =
        result.completion_coefficient * result.completion_coefficient *
        result.hyperpalinstrophy;
    result.completed_total = result.completed_negative_square +
        result.completed_enstrophy_remainder +
        result.completed_hyperpalinstrophy_remainder;
    result.completion_relative_error = relative_error(
        result.expanded_total, result.completed_total);
    const SpectralReal target_scale =
        std::pow(result.enstrophy, 1.25L) *
        std::pow(result.palinstrophy, 0.75L);
    if (target_scale > 0.0L) {
        result.absolute_target_scale_ratio =
            std::abs(result.expanded_total) / target_scale;
    }
    result.identity_verified =
        std::isfinite(result.completed_total) &&
        result.completion_relative_error < 1.0e-14L;
    return result;
}

}  // namespace lemma
