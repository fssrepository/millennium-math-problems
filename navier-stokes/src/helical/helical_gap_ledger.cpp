#include "helical_gap_ledger.hpp"

#include "triad_ledger.hpp"
#include "triad_partition.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace lemma {
namespace {

SpectralReal sector_stretching(
    const SpectralState& state,
    const std::array<SpectralIncrement, 2>& split,
    InteractionIndex interaction, std::size_t p_sign,
    std::size_t q_sign, std::size_t k_sign) {
    const auto [p_index, q_index, k_index] = interaction;
    const SpectralComplex coefficient = SpectralComplex{0.0L, 1.0L} *
        wave_dot(state.waves[q_index], split[p_sign][p_index]);
    ComplexVector advected{};
    for (std::size_t component = 0; component < 3; ++component) {
        advected[component] =
            coefficient * split[q_sign][q_index][component];
    }
    return static_cast<SpectralReal>(norm_squared(state.waves[k_index])) *
        std::real(dot_hermitian(split[k_sign][k_index], advected));
}

SpectralReal relative_residual(SpectralReal difference, SpectralReal scale) {
    return std::abs(difference) /
        std::max(std::abs(scale),
                 std::numeric_limits<SpectralReal>::min());
}

}  // namespace

HelicalGapLedgerReport HelicalGapLedger::analyze(
    const SpectralState& state) {
    HelicalGapLedgerReport report;
    std::array<SpectralIncrement, 2> split{
        SpectralIncrement(state.waves.size()),
        SpectralIncrement(state.waves.size())};
    for (std::size_t mode = 0; mode < state.waves.size(); ++mode) {
        split[0][mode] = HelicalTriadLedger::project_vector(
            state.waves[mode], state.velocity[mode], -1);
        split[1][mode] = HelicalTriadLedger::project_vector(
            state.waves[mode], state.velocity[mode], 1);
    }
    for (const InteractionIndex interaction :
         SpectralStateOps::interactions(state)) {
        const auto [p_index, q_index, k_index] = interaction;
        const int gap = TriadPartitioner::dyadic_gap(
            state.waves[p_index], state.waves[q_index],
            state.waves[k_index]);
        if (report.gaps.size() <= static_cast<std::size_t>(gap)) {
            const std::size_t previous = report.gaps.size();
            report.gaps.resize(static_cast<std::size_t>(gap + 1));
            for (std::size_t index = previous;
                 index < report.gaps.size(); ++index) {
                report.gaps[index].dyadic_gap = static_cast<int>(index);
            }
        }
        HelicalGapLedgerRow& row =
            report.gaps[static_cast<std::size_t>(gap)];
        ++row.interactions;
        for (std::size_t p_sign = 0; p_sign < 2; ++p_sign) {
            for (std::size_t q_sign = 0; q_sign < 2; ++q_sign) {
                for (std::size_t k_sign = 0; k_sign < 2; ++k_sign) {
                    const SpectralReal stretching = sector_stretching(
                        state, split, interaction,
                        p_sign, q_sign, k_sign);
                    const bool homochiral =
                        p_sign == q_sign && q_sign == k_sign;
                    if (homochiral) {
                        row.homochiral_signed += stretching;
                        row.homochiral_absolute += std::abs(stretching);
                    } else {
                        row.heterochiral_signed += stretching;
                        row.heterochiral_absolute += std::abs(stretching);
                    }
                }
            }
        }
    }

    report.homochiral_signed_total = 0.0L;
    report.heterochiral_signed_total = 0.0L;
    for (const HelicalGapLedgerRow& row : report.gaps) {
        report.homochiral_signed_total += row.homochiral_signed;
        report.heterochiral_signed_total += row.heterochiral_signed;
    }
    report.signed_total = report.homochiral_signed_total +
        report.heterochiral_signed_total;

    const TriadLedgerReport ordinary = TriadLedger::analyze(state);
    const std::size_t shared = std::min(
        report.gaps.size(), ordinary.gaps.size());
    for (std::size_t gap = 0; gap < shared; ++gap) {
        const HelicalGapLedgerRow& helical = report.gaps[gap];
        const TriadGapLedgerRow& direct = ordinary.gaps[gap];
        const SpectralReal reconstructed =
            helical.homochiral_signed + helical.heterochiral_signed;
        report.maximum_gap_reconstruction_residual = std::max(
            report.maximum_gap_reconstruction_residual,
            relative_residual(
                reconstructed - direct.signed_stretching,
                direct.absolute_pair_stretching));
    }
    report.total_reconstruction_residual = relative_residual(
        report.signed_total - ordinary.signed_total,
        ordinary.absolute_pair_total);
    return report;
}

}  // namespace lemma
