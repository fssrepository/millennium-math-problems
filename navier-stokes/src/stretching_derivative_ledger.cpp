#include "stretching_derivative_ledger.hpp"

#include <algorithm>
#include <cstddef>
#include <stdexcept>

#ifdef NS_HAVE_OPENMP
#include <omp.h>
#endif

namespace lemma {

StretchingDerivativeRoles StretchingDerivativeLedger::evaluate(
    const SpectralDynamics& dynamics,
    const SpectralState& state,
    const SpectralIncrement& direction,
    TriadSelection selection,
    int threads) {
    if (direction.size() != state.velocity.size()) {
        throw std::invalid_argument(
            "stretching derivative direction does not match state");
    }
    if (threads < 1) {
        throw std::invalid_argument(
            "stretching derivative threads must be positive");
    }
    const SpectralIncrement local_advection =
        selection.is_all()
            ? dynamics.advection(state)
            : dynamics.advection_direct_partition(state, selection);
    StretchingDerivativeRoles result;
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const SpectralReal wave2 = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        result.outer_state += wave2 * std::real(dot_hermitian(
            local_advection[mode], direction[mode]));
    }

    const std::vector<InteractionIndex>& interactions =
        SpectralStateOps::interactions(state);
    SpectralReal advecting = 0.0L;
    SpectralReal advected = 0.0L;
    const SpectralComplex imaginary_unit{0.0L, 1.0L};
#ifdef NS_HAVE_OPENMP
    const int worker_count = std::max(
        1, std::min(threads, static_cast<int>(interactions.size())));
#pragma omp parallel for num_threads(worker_count) schedule(static) \
    reduction(+ : advecting, advected)
#endif
    for (std::ptrdiff_t interaction_index = 0;
         interaction_index <
             static_cast<std::ptrdiff_t>(interactions.size());
         ++interaction_index) {
        const InteractionIndex interaction = interactions[
            static_cast<std::size_t>(interaction_index)];
        if (!TriadPartitioner::includes(state, interaction, selection)) {
            continue;
        }
        const auto [p_index, q_index, target_index] = interaction;
        const SpectralReal target_wave2 = static_cast<SpectralReal>(
            norm_squared(state.waves[target_index]));
        ComplexVector weighted_target{};
        for (std::size_t component = 0; component < 3; ++component) {
            weighted_target[component] = target_wave2 *
                state.velocity[target_index][component];
        }
        const SpectralComplex base_coefficient = imaginary_unit *
            wave_dot(state.waves[q_index], state.velocity[p_index]);
        const SpectralComplex direction_coefficient = imaginary_unit *
            wave_dot(state.waves[q_index], direction[p_index]);
        ComplexVector advecting_variation{};
        ComplexVector advected_variation{};
        for (std::size_t component = 0; component < 3; ++component) {
            advecting_variation[component] = direction_coefficient *
                state.velocity[q_index][component];
            advected_variation[component] = base_coefficient *
                direction[q_index][component];
        }
        advecting += std::real(dot_hermitian(
            weighted_target, advecting_variation));
        advected += std::real(dot_hermitian(
            weighted_target, advected_variation));
    }
    result.advecting_slot = advecting;
    result.advected_slot = advected;
    result.total = result.outer_state + result.advecting_slot +
                   result.advected_slot;
    return result;
}

}  // namespace lemma
