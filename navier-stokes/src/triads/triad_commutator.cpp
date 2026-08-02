#include "triad_commutator.hpp"

#include <algorithm>
#include <cmath>

namespace lemma {
namespace {

SpectralReal interaction_transfer(
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
    return std::real(dot_hermitian(
        state.velocity[target], contribution));
}

struct ResidualAccumulator {
    SpectralReal unweighted_residual = 0.0L;
    SpectralReal unweighted_scale = 0.0L;
    SpectralReal weighted_residual = 0.0L;
    SpectralReal weighted_scale = 0.0L;
};

void finalize_residuals(
    TriadCommutatorGapRow& row,
    const ResidualAccumulator& residual) {
    if (residual.unweighted_scale > 0.0L) {
        row.relative_unweighted_cancellation_residual =
            residual.unweighted_residual / residual.unweighted_scale;
    }
    if (residual.weighted_scale > 0.0L) {
        row.relative_weighted_identity_residual =
            residual.weighted_residual / residual.weighted_scale;
    }
}

}  // namespace

TriadCommutatorReport TriadCommutator::analyze(
    const SpectralState& state) {
    TriadCommutatorReport report;
    std::vector<ResidualAccumulator> gap_residuals;
    ResidualAccumulator total_residual;
    for (const InteractionIndex interaction :
         SpectralStateOps::interactions(state)) {
        const auto [p_index, q_index, k_index] = interaction;
        const WaveVector p = state.waves[p_index];
        const WaveVector q = state.waves[q_index];
        const WaveVector k = state.waves[k_index];
        const SpectralInteger p2 = norm_squared(p);
        const SpectralInteger q2 = norm_squared(q);
        const SpectralInteger k2 = norm_squared(k);
        if (!(p2 < q2 && p2 < k2)) {
            continue;
        }
        const WaveVector partner_advected = -k;
        if (!(q < partner_advected)) {
            continue;
        }
        const std::size_t partner_q_index =
            state.index.at(partner_advected);
        const std::size_t partner_target_index = state.index.at(-q);
        const SpectralReal first = interaction_transfer(
            state, p_index, q_index, k_index);
        const SpectralReal partner = interaction_transfer(
            state, p_index, partner_q_index, partner_target_index);
        const SpectralReal weighted_pair =
            static_cast<SpectralReal>(k2) * first +
            static_cast<SpectralReal>(q2) * partner;
        const SpectralReal commutator =
            static_cast<SpectralReal>(k2 - q2) * first;
        const SpectralReal absolute_unpaired =
            static_cast<SpectralReal>(k2) * std::abs(first) +
            static_cast<SpectralReal>(q2) * std::abs(partner);
        const int gap = TriadPartitioner::dyadic_gap(p, q, k);
        if (report.gaps.size() <= static_cast<std::size_t>(gap)) {
            const std::size_t old_size = report.gaps.size();
            report.gaps.resize(static_cast<std::size_t>(gap + 1));
            gap_residuals.resize(static_cast<std::size_t>(gap + 1));
            for (std::size_t index = old_size;
                 index < report.gaps.size(); ++index) {
                report.gaps[index].dyadic_gap =
                    static_cast<int>(index);
            }
        }
        TriadCommutatorGapRow& row =
            report.gaps[static_cast<std::size_t>(gap)];
        ResidualAccumulator& residual =
            gap_residuals[static_cast<std::size_t>(gap)];
        const SpectralReal unweighted_scale =
            std::abs(first) + std::abs(partner);
        const SpectralReal weighted_scale =
            absolute_unpaired + std::abs(commutator);
        const SpectralReal frequency_denominator =
            std::sqrt(static_cast<SpectralReal>(p2)) *
            (std::sqrt(static_cast<SpectralReal>(q2)) +
             std::sqrt(static_cast<SpectralReal>(k2)));
        const SpectralReal frequency_gain = frequency_denominator > 0.0L
            ? std::abs(static_cast<SpectralReal>(k2 - q2)) /
                  frequency_denominator
            : 0.0L;

        ++row.pairs;
        row.signed_paired_stretching += weighted_pair;
        row.absolute_unpaired_stretching += absolute_unpaired;
        row.absolute_paired_stretching += std::abs(weighted_pair);
        row.maximum_frequency_gain_ratio = std::max(
            row.maximum_frequency_gain_ratio, frequency_gain);
        residual.unweighted_residual += std::abs(first + partner);
        residual.unweighted_scale += unweighted_scale;
        residual.weighted_residual +=
            std::abs(weighted_pair - commutator);
        residual.weighted_scale += weighted_scale;

        ++report.pairs;
        report.signed_paired_stretching += weighted_pair;
        report.absolute_unpaired_stretching += absolute_unpaired;
        report.absolute_paired_stretching += std::abs(weighted_pair);
        report.maximum_frequency_gain_ratio = std::max(
            report.maximum_frequency_gain_ratio, frequency_gain);
        total_residual.unweighted_residual += std::abs(first + partner);
        total_residual.unweighted_scale += unweighted_scale;
        total_residual.weighted_residual +=
            std::abs(weighted_pair - commutator);
        total_residual.weighted_scale += weighted_scale;
    }
    for (std::size_t gap = 0; gap < report.gaps.size(); ++gap) {
        finalize_residuals(report.gaps[gap], gap_residuals[gap]);
    }
    TriadCommutatorGapRow total;
    finalize_residuals(total, total_residual);
    report.relative_unweighted_cancellation_residual =
        total.relative_unweighted_cancellation_residual;
    report.relative_weighted_identity_residual =
        total.relative_weighted_identity_residual;
    return report;
}

}  // namespace lemma
