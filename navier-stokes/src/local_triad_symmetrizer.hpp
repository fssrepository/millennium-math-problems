#pragma once

#include "spectral_state.hpp"

#include <array>
#include <cstddef>

namespace lemma {

constexpr std::size_t local_spread_bin_count = 4;

struct LocalTriadSpreadBin {
    const char* label = "";
    std::size_t triads = 0;
    SpectralReal signed_enstrophy_transfer = 0.0L;
    SpectralReal absolute_group_enstrophy_transfer = 0.0L;
    SpectralReal raw_absolute_enstrophy_transfer = 0.0L;
    SpectralReal frequency_spread_envelope = 0.0L;
};

struct LocalTriadSymmetryReport {
    std::size_t local_triads = 0;
    SpectralReal signed_local_enstrophy_transfer = 0.0L;
    SpectralReal maximum_energy_cancellation_residual = 0.0L;
    SpectralReal local_reconstruction_residual = 0.0L;
    SpectralReal maximum_frequency_spread_bound_ratio = 0.0L;
    std::array<LocalTriadSpreadBin, local_spread_bin_count> spread_bins{};
};

class LocalTriadSymmetrizer {
public:
    [[nodiscard]] static LocalTriadSymmetryReport analyze(
        const SpectralState& state);
};

}  // namespace lemma
