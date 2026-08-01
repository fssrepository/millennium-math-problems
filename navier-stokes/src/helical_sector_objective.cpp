#include "helical_sector_objective.hpp"

#include "triad_partition.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace lemma {
namespace {

std::array<SpectralIncrement, 2> split_state(const SpectralState& state) {
    std::array<SpectralIncrement, 2> split{
        SpectralIncrement(state.waves.size()),
        SpectralIncrement(state.waves.size())};
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        split[0][mode] = HelicalTriadLedger::project_vector(
            state.waves[mode], state.velocity[mode], -1);
        split[1][mode] = HelicalTriadLedger::project_vector(
            state.waves[mode], state.velocity[mode], 1);
    }
    return split;
}

void add_scaled(ComplexVector& target, const ComplexVector& source,
                SpectralComplex scale) {
    for (std::size_t component = 0; component < 3; ++component) {
        target[component] += scale * source[component];
    }
}

}  // namespace

bool HelicalSectorSelection::includes_spread(
    WaveVector first, WaveVector second, WaveVector third) const {
    const SpectralInteger minimum = std::min({
        norm_squared(first), norm_squared(second), norm_squared(third)});
    const SpectralInteger maximum = std::max({
        norm_squared(first), norm_squared(second), norm_squared(third)});
    const SpectralInteger difference = maximum - minimum;
    switch (spread) {
        case HelicalLocalSpread::all:
            return true;
        case HelicalLocalSpread::equal:
            return difference == 0;
        case HelicalLocalSpread::narrow:
            return 4 * difference <= maximum;
        case HelicalLocalSpread::broad:
            return 4 * difference > maximum;
    }
    return false;
}

HelicalSectorObjectiveValue HelicalSectorObjective::evaluate(
    const SpectralState& state, HelicalSectorSelection selection) {
    if (selection.sector_mask == 0U) {
        throw std::invalid_argument("empty helical sector selection");
    }
    HelicalSectorObjectiveValue result;
    const auto split = split_state(state);
    const SpectralComplex imaginary_unit{0.0L, 1.0L};
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const SpectralReal wave2 = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        const SpectralReal amplitude2 = std::real(dot_hermitian(
            state.velocity[mode], state.velocity[mode]));
        result.enstrophy += wave2 * amplitude2;
        result.palinstrophy += wave2 * wave2 * amplitude2;
    }
    for (const InteractionIndex interaction :
         SpectralStateOps::interactions(state)) {
        const auto [p_index, q_index, k_index] = interaction;
        if (!TriadPartitioner::is_local(
                state.waves[p_index], state.waves[q_index],
                state.waves[k_index]) ||
            !selection.includes_spread(
                state.waves[p_index], state.waves[q_index],
                state.waves[k_index])) {
            continue;
        }
        const SpectralReal k2 = static_cast<SpectralReal>(
            norm_squared(state.waves[k_index]));
        for (std::size_t p_sign = 0; p_sign < 2; ++p_sign) {
            for (std::size_t q_sign = 0; q_sign < 2; ++q_sign) {
                for (std::size_t k_sign = 0; k_sign < 2; ++k_sign) {
                    const std::size_t sector =
                        (p_sign << 2U) | (q_sign << 1U) | k_sign;
                    if (!selection.includes(sector)) {
                        continue;
                    }
                    const SpectralComplex coefficient = imaginary_unit *
                        wave_dot(state.waves[q_index],
                                 split[p_sign][p_index]);
                    ComplexVector contribution{};
                    for (std::size_t component = 0; component < 3;
                         ++component) {
                        contribution[component] = coefficient *
                            split[q_sign][q_index][component];
                    }
                    result.signed_local_stretching += k2 * std::real(
                        dot_hermitian(
                            split[k_sign][k_index], contribution));
                }
            }
        }
    }
    if (result.enstrophy > 0.0L && result.palinstrophy > 0.0L) {
        const SpectralReal stretching2 =
            result.signed_local_stretching *
            result.signed_local_stretching;
        result.critical_integrand = stretching2 * stretching2 /
            (result.enstrophy *
             result.palinstrophy * result.palinstrophy *
             result.palinstrophy);
    }
    return result;
}

SpectralIncrement HelicalSectorObjective::signed_stretching_gradient(
    const SpectralState& state, HelicalSectorSelection selection) {
    if (selection.sector_mask == 0U) {
        throw std::invalid_argument("empty helical sector selection");
    }
    const auto split = split_state(state);
    std::array<SpectralIncrement, 2> component_gradient{
        SpectralIncrement(state.waves.size()),
        SpectralIncrement(state.waves.size())};
    const SpectralComplex imaginary_unit{0.0L, 1.0L};
    const SpectralComplex minus_imaginary_unit{0.0L, -1.0L};
    for (const InteractionIndex interaction :
         SpectralStateOps::interactions(state)) {
        const auto [p_index, q_index, k_index] = interaction;
        if (!TriadPartitioner::is_local(
                state.waves[p_index], state.waves[q_index],
                state.waves[k_index]) ||
            !selection.includes_spread(
                state.waves[p_index], state.waves[q_index],
                state.waves[k_index])) {
            continue;
        }
        const SpectralReal k2 = static_cast<SpectralReal>(
            norm_squared(state.waves[k_index]));
        for (std::size_t p_sign = 0; p_sign < 2; ++p_sign) {
            for (std::size_t q_sign = 0; q_sign < 2; ++q_sign) {
                for (std::size_t k_sign = 0; k_sign < 2; ++k_sign) {
                    const std::size_t sector =
                        (p_sign << 2U) | (q_sign << 1U) | k_sign;
                    if (!selection.includes(sector)) {
                        continue;
                    }
                    const ComplexVector& advecting =
                        split[p_sign][p_index];
                    const ComplexVector& advected =
                        split[q_sign][q_index];
                    const ComplexVector& target =
                        split[k_sign][k_index];
                    const SpectralComplex coefficient = imaginary_unit *
                        wave_dot(state.waves[q_index], advecting);
                    add_scaled(
                        component_gradient[k_sign][k_index], advected,
                        k2 * coefficient);
                    add_scaled(
                        component_gradient[q_sign][q_index], target,
                        k2 * std::conj(coefficient));
                    const SpectralComplex advecting_coefficient =
                        k2 * minus_imaginary_unit *
                        dot_hermitian(advected, target);
                    const WaveVector q = state.waves[q_index];
                    const std::array<SpectralReal, 3> q_components{
                        static_cast<SpectralReal>(q.x),
                        static_cast<SpectralReal>(q.y),
                        static_cast<SpectralReal>(q.z)};
                    for (std::size_t component = 0; component < 3;
                         ++component) {
                        component_gradient[p_sign][p_index][component] +=
                            q_components[component] * advecting_coefficient;
                    }
                }
            }
        }
    }
    SpectralIncrement result(state.waves.size());
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const ComplexVector negative = HelicalTriadLedger::project_vector(
            state.waves[mode], component_gradient[0][mode], -1);
        const ComplexVector positive = HelicalTriadLedger::project_vector(
            state.waves[mode], component_gradient[1][mode], 1);
        for (std::size_t component = 0; component < 3; ++component) {
            result[mode][component] =
                negative[component] + positive[component];
        }
    }
    return result;
}

SpectralIncrement HelicalSectorObjective::critical_integrand_gradient(
    const SpectralState& state, HelicalSectorSelection selection) {
    const HelicalSectorObjectiveValue objective = evaluate(state, selection);
    SpectralIncrement result(state.waves.size());
    if (!(objective.enstrophy > 0.0L) ||
        !(objective.palinstrophy > 0.0L) ||
        objective.signed_local_stretching == 0.0L) {
        return result;
    }
    const SpectralIncrement stretching_gradient =
        signed_stretching_gradient(state, selection);
    const SpectralReal stretching = objective.signed_local_stretching;
    const SpectralReal stretching2 = stretching * stretching;
    const SpectralReal stretching4 = stretching2 * stretching2;
    const SpectralReal palinstrophy3 =
        objective.palinstrophy * objective.palinstrophy *
        objective.palinstrophy;
    const SpectralReal stretching_coefficient =
        4.0L * stretching * stretching2 /
        (objective.enstrophy * palinstrophy3);
    const SpectralReal enstrophy_coefficient =
        -stretching4 /
        (objective.enstrophy * objective.enstrophy * palinstrophy3);
    const SpectralReal palinstrophy_coefficient =
        -3.0L * stretching4 /
        (objective.enstrophy * palinstrophy3 * objective.palinstrophy);
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const SpectralReal wave2 = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        const SpectralReal wave4 = wave2 * wave2;
        for (std::size_t component = 0; component < 3; ++component) {
            result[mode][component] =
                stretching_coefficient * stretching_gradient[mode][component] +
                (2.0L * enstrophy_coefficient * wave2 +
                 2.0L * palinstrophy_coefficient * wave4) *
                    state.velocity[mode][component];
        }
    }
    return result;
}

}  // namespace lemma
