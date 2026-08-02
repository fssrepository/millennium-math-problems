#include "local_signature_objective.hpp"

#include "triad_partition.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>

namespace lemma {
namespace {

using Signature = std::array<SpectralInteger, 3>;
using SignatureTransfers = std::map<Signature, SpectralReal>;

Signature signature_of(const SpectralState& state,
                       InteractionIndex interaction) {
    const auto [p_index, q_index, k_index] = interaction;
    Signature signature{norm_squared(state.waves[p_index]),
                        norm_squared(state.waves[q_index]),
                        norm_squared(state.waves[k_index])};
    std::sort(signature.begin(), signature.end());
    return signature;
}

SpectralReal interaction_transfer(
    const SpectralState& state, InteractionIndex interaction) {
    const auto [p_index, q_index, k_index] = interaction;
    const SpectralComplex coefficient = SpectralComplex{0.0L, 1.0L} *
        wave_dot(state.waves[q_index], state.velocity[p_index]);
    ComplexVector advected{};
    for (std::size_t component = 0; component < 3; ++component) {
        advected[component] =
            coefficient * state.velocity[q_index][component];
    }
    return static_cast<SpectralReal>(
               norm_squared(state.waves[k_index])) *
        std::real(dot_hermitian(state.velocity[k_index], advected));
}

SignatureTransfers collect_transfers(const SpectralState& state) {
    SignatureTransfers transfers;
    for (const InteractionIndex interaction :
         SpectralStateOps::interactions(state)) {
        const auto [p_index, q_index, k_index] = interaction;
        if (!TriadPartitioner::is_local(
                state.waves[p_index], state.waves[q_index],
                state.waves[k_index])) {
            continue;
        }
        transfers[signature_of(state, interaction)] +=
            interaction_transfer(state, interaction);
    }
    return transfers;
}

LocalSignatureObjectiveValue summarize(
    const SignatureTransfers& transfers) {
    LocalSignatureObjectiveValue result;
    result.signatures = transfers.size();
    for (const auto& [signature, transfer] : transfers) {
        static_cast<void>(signature);
        result.signed_local_transfer += transfer;
        result.squared_signature_transfer += transfer * transfer;
    }
    if (result.squared_signature_transfer > 0.0L) {
        result.signed_amplification =
            std::abs(result.signed_local_transfer) /
            std::sqrt(result.squared_signature_transfer);
    }
    return result;
}

void add_scaled(ComplexVector& target, const ComplexVector& source,
                SpectralComplex scale) {
    for (std::size_t component = 0; component < 3; ++component) {
        target[component] += scale * source[component];
    }
}

template <typename SignatureWeight>
SpectralIncrement weighted_transfer_gradient(
    const SpectralState& state, const SignatureTransfers& transfers,
    SignatureWeight signature_weight) {
    SpectralIncrement gradient(state.waves.size());
    const SpectralComplex imaginary_unit{0.0L, 1.0L};
    const SpectralComplex minus_imaginary_unit{0.0L, -1.0L};
    for (const InteractionIndex interaction :
         SpectralStateOps::interactions(state)) {
        const auto [p_index, q_index, k_index] = interaction;
        if (!TriadPartitioner::is_local(
                state.waves[p_index], state.waves[q_index],
                state.waves[k_index])) {
            continue;
        }
        const SpectralReal weight = signature_weight(
            transfers.at(signature_of(state, interaction)));
        const SpectralReal k2 = static_cast<SpectralReal>(
            norm_squared(state.waves[k_index]));
        const ComplexVector& advecting = state.velocity[p_index];
        const ComplexVector& advected = state.velocity[q_index];
        const ComplexVector& target = state.velocity[k_index];
        const SpectralComplex coefficient = imaginary_unit *
            wave_dot(state.waves[q_index], advecting);
        add_scaled(gradient[k_index], advected,
                   weight * k2 * coefficient);
        add_scaled(gradient[q_index], target,
                   weight * k2 * std::conj(coefficient));
        const SpectralComplex advecting_coefficient =
            weight * k2 * minus_imaginary_unit *
            dot_hermitian(advected, target);
        const WaveVector q = state.waves[q_index];
        const std::array<SpectralReal, 3> q_components{
            static_cast<SpectralReal>(q.x),
            static_cast<SpectralReal>(q.y),
            static_cast<SpectralReal>(q.z)};
        for (std::size_t component = 0; component < 3; ++component) {
            gradient[p_index][component] +=
                q_components[component] * advecting_coefficient;
        }
    }
    return gradient;
}

}  // namespace

LocalSignatureObjectiveValue LocalSignatureObjective::evaluate(
    const SpectralState& state) {
    return summarize(collect_transfers(state));
}

SpectralIncrement LocalSignatureObjective::signed_amplification_gradient(
    const SpectralState& state) {
    const SignatureTransfers transfers = collect_transfers(state);
    const LocalSignatureObjectiveValue objective = summarize(transfers);
    if (!(objective.squared_signature_transfer > 0.0L) ||
        objective.signed_local_transfer == 0.0L) {
        return SpectralIncrement(state.waves.size());
    }

    const SpectralReal root_square_sum =
        std::sqrt(objective.squared_signature_transfer);
    const SpectralReal inverse_root = 1.0L / root_square_sum;
    const SpectralReal denominator3 =
        objective.squared_signature_transfer * root_square_sum;
    const SpectralReal signed_coefficient =
        std::copysign(inverse_root, objective.signed_local_transfer);
    return weighted_transfer_gradient(
        state, transfers, [&](SpectralReal signature_transfer) {
            return signed_coefficient -
            std::abs(objective.signed_local_transfer) *
                signature_transfer / denominator3;
        });
}

SpectralIncrement LocalSignatureObjective::absolute_signed_transfer_gradient(
    const SpectralState& state) {
    const SignatureTransfers transfers = collect_transfers(state);
    const LocalSignatureObjectiveValue objective = summarize(transfers);
    if (objective.signed_local_transfer == 0.0L) {
        return SpectralIncrement(state.waves.size());
    }
    const SpectralReal sign = std::copysign(
        1.0L, objective.signed_local_transfer);
    return weighted_transfer_gradient(
        state, transfers, [=](SpectralReal) { return sign; });
}

}  // namespace lemma
