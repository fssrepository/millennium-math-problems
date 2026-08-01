#include "spectral_objective.hpp"

#include <cmath>
#include <cstddef>

namespace lemma {

SpectralObjective::SpectralObjective(const SpectralDynamics& dynamics)
    : dynamics_(dynamics) {}

StaticObjective SpectralObjective::evaluate(
    const SpectralState& state, TriadPartition partition) const {
    const SpectralIncrement advection = partition == TriadPartition::all
        ? dynamics_.advection(state)
        : dynamics_.advection_direct_partition(state, partition);
    StaticObjective result;
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const SpectralReal wave2 =
            static_cast<SpectralReal>(norm_squared(state.waves[mode]));
        const SpectralReal mode_energy = std::real(
            dot_hermitian(state.velocity[mode], state.velocity[mode]));
        const SpectralReal pairing = std::real(
            dot_hermitian(state.velocity[mode], advection[mode]));
        result.energy += mode_energy;
        result.enstrophy += wave2 * mode_energy;
        result.palinstrophy += wave2 * wave2 * mode_energy;
        result.signed_vortex_stretching += wave2 * pairing;
    }
    result.vortex_stretching =
        std::abs(result.signed_vortex_stretching);
    const SpectralReal denominator =
        std::pow(result.enstrophy, 0.75L) *
        std::pow(result.palinstrophy, 0.75L);
    if (denominator > 0.0L) {
        result.depletion = result.vortex_stretching / denominator;
        const SpectralReal depletion4 =
            std::pow(result.depletion, 4.0L);
        result.energy_level_quantity = depletion4 * result.enstrophy;
        result.critical_integrand =
            result.energy_level_quantity * result.enstrophy;
    }
    return result;
}

SpectralIncrement SpectralObjective::energy_level_gradient(
    const SpectralState& state) const {
    const StaticObjective objective = evaluate(state);
    SpectralIncrement result(state.waves.size());
    if (!(objective.enstrophy > 0.0L) ||
        !(objective.palinstrophy > 0.0L) ||
        objective.signed_vortex_stretching == 0.0L) {
        return result;
    }

    const SpectralIncrement stretching_gradient =
        signed_stretching_gradient(state, TriadPartition::all);

    const SpectralReal stretching = objective.signed_vortex_stretching;
    const SpectralReal stretching2 = stretching * stretching;
    const SpectralReal stretching4 = stretching2 * stretching2;
    const SpectralReal enstrophy2 =
        objective.enstrophy * objective.enstrophy;
    const SpectralReal palinstrophy2 =
        objective.palinstrophy * objective.palinstrophy;
    const SpectralReal palinstrophy3 =
        palinstrophy2 * objective.palinstrophy;
    const SpectralReal stretching_coefficient =
        4.0L * stretching * stretching2 /
        (enstrophy2 * palinstrophy3);
    const SpectralReal enstrophy_coefficient =
        -2.0L * stretching4 /
        (enstrophy2 * objective.enstrophy * palinstrophy3);
    const SpectralReal palinstrophy_coefficient =
        -3.0L * stretching4 /
        (enstrophy2 * palinstrophy3 * objective.palinstrophy);
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const SpectralReal wave2 =
            static_cast<SpectralReal>(norm_squared(state.waves[mode]));
        const SpectralReal wave4 = wave2 * wave2;
        for (std::size_t component = 0; component < 3; ++component) {
            result[mode][component] =
                stretching_coefficient *
                    stretching_gradient[mode][component] +
                (2.0L * enstrophy_coefficient * wave2 +
                 2.0L * palinstrophy_coefficient * wave4) *
                    state.velocity[mode][component];
        }
    }
    return result;
}

SpectralIncrement SpectralObjective::critical_integrand_gradient(
    const SpectralState& state, TriadPartition partition) const {
    const StaticObjective objective = evaluate(state, partition);
    SpectralIncrement result(state.waves.size());
    if (!(objective.enstrophy > 0.0L) ||
        !(objective.palinstrophy > 0.0L) ||
        objective.signed_vortex_stretching == 0.0L) {
        return result;
    }
    const SpectralIncrement stretching_gradient =
        signed_stretching_gradient(state, partition);
    const SpectralReal stretching = objective.signed_vortex_stretching;
    const SpectralReal stretching2 = stretching * stretching;
    const SpectralReal stretching4 = stretching2 * stretching2;
    const SpectralReal palinstrophy2 =
        objective.palinstrophy * objective.palinstrophy;
    const SpectralReal palinstrophy3 =
        palinstrophy2 * objective.palinstrophy;
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
        const SpectralReal wave2 =
            static_cast<SpectralReal>(norm_squared(state.waves[mode]));
        const SpectralReal wave4 = wave2 * wave2;
        for (std::size_t component = 0; component < 3; ++component) {
            result[mode][component] =
                stretching_coefficient *
                    stretching_gradient[mode][component] +
                (2.0L * enstrophy_coefficient * wave2 +
                 2.0L * palinstrophy_coefficient * wave4) *
                    state.velocity[mode][component];
        }
    }
    return result;
}

SpectralIncrement SpectralObjective::signed_stretching_gradient(
    const SpectralState& state, TriadPartition partition) const {
    const SpectralIncrement advection = partition == TriadPartition::all
        ? dynamics_.advection(state)
        : dynamics_.advection_direct_partition(state, partition);
    SpectralIncrement weighted_velocity(state.waves.size());
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const SpectralReal wave2 =
            static_cast<SpectralReal>(norm_squared(state.waves[mode]));
        for (std::size_t component = 0; component < 3; ++component) {
            weighted_velocity[mode][component] =
                wave2 * state.velocity[mode][component];
        }
    }
    SpectralIncrement result = partition == TriadPartition::all
        ? dynamics_.advection_vjp(state, weighted_velocity)
        : dynamics_.advection_vjp_direct_partition(
              state, weighted_velocity, partition);
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const SpectralReal wave2 =
            static_cast<SpectralReal>(norm_squared(state.waves[mode]));
        for (std::size_t component = 0; component < 3; ++component) {
            result[mode][component] +=
                wave2 * advection[mode][component];
        }
    }
    return result;
}

}  // namespace lemma
