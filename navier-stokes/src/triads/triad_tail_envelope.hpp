#pragma once

#include "triad_ledger.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace lemma {

constexpr std::size_t separated_low_role_count = 3;

struct TriadTailEnvelopeGapRow {
    int dyadic_gap = 0;
    std::array<std::size_t, separated_low_role_count> terms_by_low_role{};
    std::array<SpectralReal, separated_low_role_count>
        signed_stretching_by_low_role{};
    std::array<SpectralReal, separated_low_role_count>
        absolute_stretching_by_low_role{};
    std::array<SpectralReal, separated_low_role_count>
        amplitude_envelope_by_low_role{};
    std::array<SpectralReal, separated_low_role_count>
        maximum_amplitude_bound_ratio{};
    std::array<SpectralReal, separated_low_role_count>
        maximum_normalized_frequency_ratio{};
};

struct TriadTailEnvelopeReport {
    std::array<std::size_t, separated_low_role_count> terms_by_low_role{};
    std::array<SpectralReal, separated_low_role_count>
        signed_stretching_by_low_role{};
    std::array<SpectralReal, separated_low_role_count>
        absolute_stretching_by_low_role{};
    std::array<SpectralReal, separated_low_role_count>
        amplitude_envelope_by_low_role{};
    std::array<SpectralReal, separated_low_role_count>
        maximum_amplitude_bound_ratio{};
    std::array<SpectralReal, separated_low_role_count>
        maximum_normalized_frequency_ratio{};
    std::vector<TriadTailEnvelopeGapRow> gaps;
};

class TriadTailEnvelope {
public:
    // For every separated triad, place the unique low wave in one role.
    // Low-advecting terms use the exact commutator pair. Low-advected terms
    // already carry a low derivative, and low-target terms carry a low
    // enstrophy weight. The normalized frequency ratios are bounded by one.
    [[nodiscard]] static TriadTailEnvelopeReport analyze(
        const SpectralState& state);
};

}  // namespace lemma
