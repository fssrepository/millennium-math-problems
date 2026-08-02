#include "triad_tail_envelope.hpp"

#include <algorithm>
#include <cmath>

namespace lemma {
namespace {

struct InteractionMeasurement {
    SpectralReal transfer = 0.0L;
    SpectralReal amplitude_product = 0.0L;
};

SpectralReal vector_norm(const ComplexVector& value) {
    return std::sqrt(std::max(
        0.0L, std::real(dot_hermitian(value, value))));
}

InteractionMeasurement measure_interaction(
    const SpectralState& state, std::size_t advecting,
    std::size_t advected, std::size_t target) {
    const SpectralComplex imaginary_unit{0.0L, 1.0L};
    const SpectralComplex coefficient = imaginary_unit * wave_dot(
        state.waves[advected], state.velocity[advecting]);
    ComplexVector contribution{};
    for (std::size_t component = 0; component < 3; ++component) {
        contribution[component] =
            coefficient * state.velocity[advected][component];
    }
    return {
        std::real(dot_hermitian(state.velocity[target], contribution)),
        vector_norm(state.velocity[advecting]) *
            vector_norm(state.velocity[advected]) *
            vector_norm(state.velocity[target])};
}

void accumulate(
    TriadTailEnvelopeGapRow& row, TriadLowRole role,
    SpectralReal signed_stretching, SpectralReal envelope,
    SpectralReal normalized_frequency_ratio) {
    const std::size_t role_index = static_cast<std::size_t>(role);
    ++row.terms_by_low_role[role_index];
    row.signed_stretching_by_low_role[role_index] += signed_stretching;
    row.absolute_stretching_by_low_role[role_index] +=
        std::abs(signed_stretching);
    row.amplitude_envelope_by_low_role[role_index] += envelope;
    if (envelope > 0.0L) {
        row.maximum_amplitude_bound_ratio[role_index] = std::max(
            row.maximum_amplitude_bound_ratio[role_index],
            std::abs(signed_stretching) / envelope);
    }
    row.maximum_normalized_frequency_ratio[role_index] = std::max(
        row.maximum_normalized_frequency_ratio[role_index],
        normalized_frequency_ratio);
}

void merge(
    TriadTailEnvelopeReport& report,
    const TriadTailEnvelopeGapRow& row) {
    for (std::size_t role = 0; role < separated_low_role_count; ++role) {
        report.terms_by_low_role[role] += row.terms_by_low_role[role];
        report.signed_stretching_by_low_role[role] +=
            row.signed_stretching_by_low_role[role];
        report.absolute_stretching_by_low_role[role] +=
            row.absolute_stretching_by_low_role[role];
        report.amplitude_envelope_by_low_role[role] +=
            row.amplitude_envelope_by_low_role[role];
        report.maximum_amplitude_bound_ratio[role] = std::max(
            report.maximum_amplitude_bound_ratio[role],
            row.maximum_amplitude_bound_ratio[role]);
        report.maximum_normalized_frequency_ratio[role] = std::max(
            report.maximum_normalized_frequency_ratio[role],
            row.maximum_normalized_frequency_ratio[role]);
    }
}

}  // namespace

TriadTailEnvelopeReport TriadTailEnvelope::analyze(
    const SpectralState& state) {
    TriadTailEnvelopeReport report;
    for (const InteractionIndex interaction :
         SpectralStateOps::interactions(state)) {
        const auto [p_index, q_index, k_index] = interaction;
        const WaveVector p = state.waves[p_index];
        const WaveVector q = state.waves[q_index];
        const WaveVector k = state.waves[k_index];
        const SpectralInteger p2 = norm_squared(p);
        const SpectralInteger q2 = norm_squared(q);
        const SpectralInteger k2 = norm_squared(k);
        const SpectralInteger low2 = std::min({p2, q2, k2});
        const int low_count = static_cast<int>(p2 == low2) +
                              static_cast<int>(q2 == low2) +
                              static_cast<int>(k2 == low2);
        if (low_count != 1) {
            continue;
        }
        const int gap = TriadPartitioner::dyadic_gap(p, q, k);
        if (gap == 0) {
            continue;
        }
        if (report.gaps.size() <= static_cast<std::size_t>(gap)) {
            const std::size_t old_size = report.gaps.size();
            report.gaps.resize(static_cast<std::size_t>(gap + 1));
            for (std::size_t index = old_size;
                 index < report.gaps.size(); ++index) {
                report.gaps[index].dyadic_gap =
                    static_cast<int>(index);
            }
        }
        TriadTailEnvelopeGapRow& row =
            report.gaps[static_cast<std::size_t>(gap)];
        const SpectralReal p_norm =
            std::sqrt(static_cast<SpectralReal>(p2));
        const SpectralReal q_norm =
            std::sqrt(static_cast<SpectralReal>(q2));
        const SpectralReal k_norm =
            std::sqrt(static_cast<SpectralReal>(k2));

        if (p2 == low2) {
            const WaveVector partner_advected = -k;
            if (!(q < partner_advected)) {
                continue;
            }
            const std::size_t partner_q_index =
                state.index.at(partner_advected);
            const std::size_t partner_target_index = state.index.at(-q);
            const InteractionMeasurement first = measure_interaction(
                state, p_index, q_index, k_index);
            const InteractionMeasurement partner = measure_interaction(
                state, p_index, partner_q_index,
                partner_target_index);
            const SpectralReal signed_stretching =
                static_cast<SpectralReal>(k2) * first.transfer +
                static_cast<SpectralReal>(q2) * partner.transfer;
            const SpectralReal frequency_weight =
                std::abs(static_cast<SpectralReal>(k2 - q2)) * q_norm;
            const SpectralReal envelope =
                frequency_weight * first.amplitude_product;
            const SpectralReal high = std::max(q_norm, k_norm);
            const SpectralReal normalized_frequency_ratio =
                high > 0.0L && p_norm > 0.0L
                    ? frequency_weight /
                          (2.0L * p_norm * high * high)
                    : 0.0L;
            accumulate(row, TriadLowRole::advecting,
                       signed_stretching, envelope,
                       normalized_frequency_ratio);
            static_cast<void>(partner.amplitude_product);
            continue;
        }

        const InteractionMeasurement measurement = measure_interaction(
            state, p_index, q_index, k_index);
        const SpectralReal signed_stretching =
            static_cast<SpectralReal>(k2) * measurement.transfer;
        if (q2 == low2) {
            const SpectralReal frequency_weight =
                static_cast<SpectralReal>(k2) * q_norm;
            const SpectralReal high = std::max(p_norm, k_norm);
            const SpectralReal normalized_frequency_ratio =
                high > 0.0L && q_norm > 0.0L
                    ? frequency_weight / (q_norm * high * high)
                    : 0.0L;
            accumulate(
                row, TriadLowRole::advected, signed_stretching,
                frequency_weight * measurement.amplitude_product,
                normalized_frequency_ratio);
            continue;
        }

        const SpectralReal frequency_weight =
            static_cast<SpectralReal>(k2) * q_norm;
        const SpectralReal high = std::max(p_norm, q_norm);
        const SpectralReal normalized_frequency_ratio =
            high > 0.0L && k_norm > 0.0L
                ? frequency_weight /
                      (k_norm * k_norm * high)
                : 0.0L;
        accumulate(
            row, TriadLowRole::target, signed_stretching,
            frequency_weight * measurement.amplitude_product,
            normalized_frequency_ratio);
    }
    for (const TriadTailEnvelopeGapRow& row : report.gaps) {
        merge(report, row);
    }
    return report;
}

}  // namespace lemma
