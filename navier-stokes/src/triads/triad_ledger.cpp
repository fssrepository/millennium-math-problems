#include "triad_ledger.hpp"

#include <algorithm>
#include <cmath>

namespace lemma {
namespace {

TriadLowRole classify_low_role(
    SpectralInteger advecting, SpectralInteger advected,
    SpectralInteger target) {
    const SpectralInteger smallest =
        std::min({advecting, advected, target});
    const int minima = static_cast<int>(advecting == smallest) +
                       static_cast<int>(advected == smallest) +
                       static_cast<int>(target == smallest);
    if (minima != 1) {
        return TriadLowRole::tied;
    }
    if (advecting == smallest) {
        return TriadLowRole::advecting;
    }
    if (advected == smallest) {
        return TriadLowRole::advected;
    }
    return TriadLowRole::target;
}

}  // namespace

TriadLedgerReport TriadLedger::analyze(const SpectralState& state) {
    TriadLedgerReport report;
    const SpectralComplex imaginary_unit{0.0L, 1.0L};
    for (const InteractionIndex interaction :
         SpectralStateOps::interactions(state)) {
        const auto [p_index, q_index, target_index] = interaction;
        const WaveVector p = state.waves[p_index];
        const WaveVector q = state.waves[q_index];
        const WaveVector k = state.waves[target_index];
        const int gap = TriadPartitioner::dyadic_gap(p, q, k);
        if (report.gaps.size() <= static_cast<std::size_t>(gap)) {
            const std::size_t old_size = report.gaps.size();
            report.gaps.resize(static_cast<std::size_t>(gap + 1));
            for (std::size_t index = old_size;
                 index < report.gaps.size(); ++index) {
                report.gaps[index].dyadic_gap =
                    static_cast<int>(index);
            }
        }
        TriadGapLedgerRow& row =
            report.gaps[static_cast<std::size_t>(gap)];
        const SpectralComplex coefficient =
            imaginary_unit * wave_dot(q, state.velocity[p_index]);
        ComplexVector pair{};
        for (std::size_t component = 0; component < 3; ++component) {
            pair[component] =
                coefficient * state.velocity[q_index][component];
        }
        const SpectralReal stretching =
            static_cast<SpectralReal>(norm_squared(k)) *
            std::real(dot_hermitian(
                state.velocity[target_index], pair));
        const SpectralReal absolute_stretching = std::abs(stretching);
        const TriadLowRole role = classify_low_role(
            norm_squared(p), norm_squared(q), norm_squared(k));
        const std::size_t role_index = static_cast<std::size_t>(role);

        ++row.interactions;
        row.signed_stretching += stretching;
        row.absolute_pair_stretching += absolute_stretching;
        ++row.interactions_by_low_role[role_index];
        row.signed_stretching_by_low_role[role_index] += stretching;
        row.absolute_stretching_by_low_role[role_index] +=
            absolute_stretching;
        report.signed_total += stretching;
        report.absolute_pair_total += absolute_stretching;
        if (gap == 0) {
            report.signed_local += stretching;
        } else {
            report.signed_nonlocal += stretching;
        }
    }
    report.partition_residual = std::abs(
        report.signed_total - report.signed_local -
        report.signed_nonlocal);
    return report;
}

const char* TriadLedger::low_role_name(TriadLowRole role) {
    switch (role) {
        case TriadLowRole::advecting:
            return "advecting";
        case TriadLowRole::advected:
            return "advected";
        case TriadLowRole::target:
            return "target";
        case TriadLowRole::tied:
            return "tied";
    }
    return "unknown";
}

}  // namespace lemma
