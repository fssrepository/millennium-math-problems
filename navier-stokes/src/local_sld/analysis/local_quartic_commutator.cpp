#include "local_quartic_commutator.hpp"

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

LocalQuarticCommutatorReport LocalQuarticCommutator::evaluate(
    const SpectralDynamics& dynamics,
    const SpectralState& state,
    const LocalCriticalDerivativeLedgerReport& derivative_ledger) {
    const SpectralIncrement local_advection =
        dynamics.advection_direct_partition(
            state, TriadPartition::local);
    SpectralIncrement weighted_velocity(state.waves.size());
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const SpectralReal wave2 = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        for (std::size_t component = 0; component < 3; ++component) {
            weighted_velocity[mode][component] =
                wave2 * state.velocity[mode][component];
        }
    }
    const SpectralIncrement transported_weighted_velocity =
        dynamics.advection_bilinear_direct_partition(
            state, state.velocity, weighted_velocity,
            TriadPartition::local);

    LocalQuarticCommutatorReport result;
    result.outer_state_derivative = derivative_ledger
        .stretching_local_nonlinear_roles.outer_state;
    result.advected_slot_derivative = derivative_ledger
        .stretching_local_nonlinear_roles.advected_slot;
    result.combined_derivative = result.outer_state_derivative +
        result.advected_slot_derivative;
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        const SpectralReal wave2 = static_cast<SpectralReal>(
            norm_squared(state.waves[mode]));
        ComplexVector commutator{};
        for (std::size_t component = 0; component < 3; ++component) {
            commutator[component] =
                wave2 * local_advection[mode][component] -
                transported_weighted_velocity[mode][component];
        }
        result.negative_commutator_pairing -= std::real(
            dot_hermitian(local_advection[mode], commutator));
        result.local_advection_l2_squared += std::real(dot_hermitian(
            local_advection[mode], local_advection[mode]));
        result.commutator_l2_squared += std::real(
            dot_hermitian(commutator, commutator));
    }
    const SpectralReal cauchy_denominator = std::sqrt(
        result.local_advection_l2_squared *
        result.commutator_l2_squared);
    if (cauchy_denominator > 0.0L) {
        result.cauchy_ratio =
            std::abs(result.negative_commutator_pairing) /
            cauchy_denominator;
    }
    result.identity_relative_error = relative_error(
        result.combined_derivative,
        result.negative_commutator_pairing);

    for (const InteractionIndex interaction :
         SpectralStateOps::interactions(state)) {
        if (!TriadPartitioner::includes(
                state, interaction, TriadPartition::local)) {
            continue;
        }
        const auto [p_index, q_index, target_index] = interaction;
        const WaveVector p = state.waves[p_index];
        const WaveVector q = state.waves[q_index];
        const WaveVector k = state.waves[target_index];
        const SpectralReal numerator = std::abs(
            static_cast<SpectralReal>(
                norm_squared(k) - norm_squared(q)));
        const SpectralReal denominator = std::sqrt(
            static_cast<SpectralReal>(norm_squared(p))) *
            (std::sqrt(static_cast<SpectralReal>(norm_squared(k))) +
             std::sqrt(static_cast<SpectralReal>(norm_squared(q))));
        if (denominator > 0.0L) {
            result.maximum_symbol_ratio = std::max(
                result.maximum_symbol_ratio,
                numerator / denominator);
        }
    }
    // Equality is attained by collinear integer wave vectors.  The square
    // roots in the diagnostic ratio can round a few ulps above one even
    // though |k^2-q^2|=|p.(k+q)| <= |p|(|k|+|q|) holds exactly.
    constexpr SpectralReal floating_point_slack = 1.0e-15L;
    result.symbol_bound_holds =
        result.maximum_symbol_ratio <= 1.0L + floating_point_slack;
    result.finite =
        std::isfinite(result.negative_commutator_pairing) &&
        result.identity_relative_error < 1.0e-12L &&
        result.cauchy_ratio <= 1.0L + floating_point_slack &&
        result.symbol_bound_holds;
    return result;
}

}  // namespace lemma
