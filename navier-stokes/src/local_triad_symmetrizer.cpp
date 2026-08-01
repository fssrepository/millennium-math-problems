#include "local_triad_symmetrizer.hpp"

#include "triad_ledger.hpp"
#include "triad_partition.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>

namespace lemma {
namespace {

using TriadKey = std::array<WaveVector, 3>;

struct TriadAccumulator {
    std::array<SpectralReal, 3> target_energy_transfer{};
    std::array<SpectralReal, 3> target_absolute_transfer{};
};

TriadKey interaction_key(
    WaveVector p, WaveVector q, WaveVector k) {
    TriadKey key{-k, p, q};
    std::sort(key.begin(), key.end());
    return key;
}

std::size_t target_slot(const TriadKey& key, WaveVector target) {
    const auto found = std::find(key.begin(), key.end(), -target);
    return static_cast<std::size_t>(found - key.begin());
}

SpectralReal energy_transfer(
    const SpectralState& state, InteractionIndex interaction) {
    const auto [p_index, q_index, k_index] = interaction;
    const SpectralComplex coefficient = SpectralComplex{0.0L, 1.0L} *
        wave_dot(state.waves[q_index], state.velocity[p_index]);
    ComplexVector advected{};
    for (std::size_t component = 0; component < 3; ++component) {
        advected[component] =
            coefficient * state.velocity[q_index][component];
    }
    return std::real(dot_hermitian(
        state.velocity[k_index], advected));
}

std::size_t spread_bin(
    SpectralInteger minimum_wave2, SpectralInteger maximum_wave2) {
    const SpectralInteger difference = maximum_wave2 - minimum_wave2;
    if (difference == 0) {
        return 0;
    }
    if (4 * difference <= maximum_wave2) {
        return 1;
    }
    if (2 * difference <= maximum_wave2) {
        return 2;
    }
    return 3;
}

SpectralReal relative_residual(SpectralReal difference, SpectralReal scale) {
    return std::abs(difference) /
        std::max(std::abs(scale),
                 std::numeric_limits<SpectralReal>::min());
}

}  // namespace

LocalTriadSymmetryReport LocalTriadSymmetrizer::analyze(
    const SpectralState& state) {
    std::map<TriadKey, TriadAccumulator> groups;
    for (const InteractionIndex interaction :
         SpectralStateOps::interactions(state)) {
        const auto [p_index, q_index, k_index] = interaction;
        const WaveVector p = state.waves[p_index];
        const WaveVector q = state.waves[q_index];
        const WaveVector k = state.waves[k_index];
        if (!TriadPartitioner::is_local(p, q, k)) {
            continue;
        }
        const TriadKey key = interaction_key(p, q, k);
        const std::size_t slot = target_slot(key, k);
        if (slot >= key.size()) {
            continue;
        }
        const SpectralReal transfer = energy_transfer(state, interaction);
        groups[key].target_energy_transfer[slot] += transfer;
        groups[key].target_absolute_transfer[slot] += std::abs(transfer);
    }

    LocalTriadSymmetryReport report;
    report.spread_bins[0].label = "equal";
    report.spread_bins[1].label = "relative-spread<=1/4";
    report.spread_bins[2].label = "relative-spread<=1/2";
    report.spread_bins[3].label = "relative-spread>1/2";
    SpectralReal total_absolute_energy_transfer = 0.0L;
    for (const auto& [key, group] : groups) {
        static_cast<void>(key);
        for (const SpectralReal value : group.target_absolute_transfer) {
            total_absolute_energy_transfer += value;
        }
    }
    SpectralReal reconstruction_scale = 0.0L;
    for (const auto& [key, group] : groups) {
        const std::array<SpectralInteger, 3> wave2{
            norm_squared(key[0]), norm_squared(key[1]),
            norm_squared(key[2])};
        const SpectralInteger minimum_wave2 = *std::min_element(
            wave2.begin(), wave2.end());
        const SpectralInteger maximum_wave2 = *std::max_element(
            wave2.begin(), wave2.end());
        SpectralReal energy_sum = 0.0L;
        SpectralReal energy_absolute = 0.0L;
        SpectralReal weighted = 0.0L;
        SpectralReal raw_weighted_absolute = 0.0L;
        for (std::size_t target = 0; target < 3; ++target) {
            energy_sum += group.target_energy_transfer[target];
            energy_absolute += group.target_absolute_transfer[target];
            weighted += static_cast<SpectralReal>(wave2[target]) *
                group.target_energy_transfer[target];
            raw_weighted_absolute +=
                static_cast<SpectralReal>(wave2[target]) *
                group.target_absolute_transfer[target];
        }
        if (energy_absolute >
            1e-14L * total_absolute_energy_transfer) {
            report.maximum_energy_cancellation_residual = std::max(
                report.maximum_energy_cancellation_residual,
                relative_residual(energy_sum, energy_absolute));
        }
        const std::size_t bin = spread_bin(
            minimum_wave2, maximum_wave2);
        const SpectralReal centered_weighted = weighted -
            static_cast<SpectralReal>(minimum_wave2) * energy_sum;
        const SpectralReal spread_envelope =
            static_cast<SpectralReal>(maximum_wave2 - minimum_wave2) *
            energy_absolute;
        LocalTriadSpreadBin& row = report.spread_bins[bin];
        ++row.triads;
        row.signed_enstrophy_transfer += weighted;
        row.absolute_group_enstrophy_transfer += std::abs(weighted);
        row.raw_absolute_enstrophy_transfer += raw_weighted_absolute;
        row.frequency_spread_envelope += spread_envelope;
        if (spread_envelope >
            1e-14L * total_absolute_energy_transfer) {
            report.maximum_frequency_spread_bound_ratio = std::max(
                report.maximum_frequency_spread_bound_ratio,
                std::abs(centered_weighted) / spread_envelope);
        }
        ++report.local_triads;
        report.signed_local_enstrophy_transfer += weighted;
        reconstruction_scale += raw_weighted_absolute;
    }

    const TriadLedgerReport direct = TriadLedger::analyze(state);
    report.local_reconstruction_residual = relative_residual(
        report.signed_local_enstrophy_transfer - direct.signed_local,
        reconstruction_scale);
    return report;
}

}  // namespace lemma
