#pragma once

#include "spectral_state.hpp"

#include <cstdint>

namespace lemma {

struct TriadCertificate {
    int modes = 0;
    int samples = 0;
    SpectralReal maximum_normalized_energy_residual = 0.0L;
    SpectralReal maximum_divergence_residual = 0.0L;
    SpectralReal maximum_reality_residual = 0.0L;
    SpectralReal maximum_classical_ratio = 0.0L;
    SpectralReal maximum_vortex_stretching = 0.0L;
    SpectralReal maximum_detailed_triad_residual = 0.0L;
    SpectralReal maximum_relative_detailed_triad_residual = 0.0L;
    SpectralReal maximum_nonlocal_absolute_fraction = 0.0L;
    SpectralReal maximum_flux_efficiency = 0.0L;
    SpectralReal maximum_local_cumulative_flux = 0.0L;
    SpectralReal maximum_nonlocal_cumulative_flux = 0.0L;
    SpectralReal maximum_flux_partition_residual = 0.0L;
    bool nonzero_vortex_stretching_seen = false;
};

class TriadVerifier {
public:
    [[nodiscard]] static TriadCertificate analyze(
        int cutoff, int samples, std::uint64_t seed);

private:
    struct InteractionAnalysis;
    struct TriadMeasurement;

    [[nodiscard]] static InteractionAnalysis analyze_interactions(
        const SpectralState& state);
    [[nodiscard]] static TriadMeasurement measure(
        const SpectralState& state);
};

}  // namespace lemma
