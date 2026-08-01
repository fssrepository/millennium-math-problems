#pragma once

#include "spectral_state.hpp"

#include <array>
#include <cstddef>

namespace lemma {

constexpr std::size_t helical_sector_count = 8;

struct HelicalSectorRow {
    int advecting_sign = -1;
    int advected_sign = -1;
    int target_sign = -1;
    SpectralReal signed_total_stretching = 0.0L;
    SpectralReal absolute_total_stretching = 0.0L;
    SpectralReal signed_local_stretching = 0.0L;
    SpectralReal absolute_local_stretching = 0.0L;
};

struct HelicalTriadReport {
    SpectralReal positive_helical_energy = 0.0L;
    SpectralReal negative_helical_energy = 0.0L;
    SpectralReal helicity = 0.0L;
    SpectralReal relative_velocity_reconstruction_residual = 0.0L;
    SpectralReal signed_total_stretching = 0.0L;
    SpectralReal signed_local_stretching = 0.0L;
    SpectralReal homochiral_local_stretching = 0.0L;
    SpectralReal heterochiral_local_stretching = 0.0L;
    SpectralReal homochiral_absolute_local_stretching = 0.0L;
    SpectralReal heterochiral_absolute_local_stretching = 0.0L;
    SpectralReal relative_total_reconstruction_residual = 0.0L;
    SpectralReal relative_local_reconstruction_residual = 0.0L;
    std::array<HelicalSectorRow, helical_sector_count> sectors{};
};

class HelicalTriadLedger {
public:
    [[nodiscard]] static HelicalTriadReport analyze(
        const SpectralState& state);
};

}  // namespace lemma
