#pragma once

#include "periodic_shell_geometry.hpp"
#include "spectral_state.hpp"

#include <cstdint>

namespace lemma {

struct PeriodicTailBoundEvaluation {
    int cutoff = 0;
    int minimum_gap = 0;
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal signed_tail_stretching = 0.0L;
    SpectralReal absolute_signed_tail_stretching = 0.0L;
    SpectralReal one_gain_bound = 0.0L;
    SpectralReal three_gain_bound = 0.0L;
    SpectralReal total_bound = 0.0L;
    SpectralReal bound_ratio = 0.0L;
};

struct PeriodicTailBoundCertificate {
    int cutoff = 0;
    int minimum_gap = 0;
    int samples = 0;
    std::uint64_t seed = 0;
    SpectralReal maximum_bound_ratio = 0.0L;
    bool nonzero_tail_seen = false;
    bool all_bounds_hold = true;
};

class PeriodicTailBound {
public:
    [[nodiscard]] static PeriodicTailBoundEvaluation evaluate(
        const SpectralState& state, int minimum_gap,
        const PeriodicShellGeometryCertificate& geometry);

    [[nodiscard]] static PeriodicTailBoundCertificate verify_random(
        int cutoff, int minimum_gap, int samples, std::uint64_t seed,
        const PeriodicShellGeometryCertificate& geometry);
};

}  // namespace lemma
