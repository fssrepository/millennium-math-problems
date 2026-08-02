#pragma once

#include "moving_gap_controller.hpp"
#include "periodic_shell_geometry.hpp"

#include <cstdint>

namespace lemma {

struct FarTailClosureEvaluation {
    SpectralReal viscosity = 0.0L;
    SpectralReal enstrophy = 0.0L;
    MovingGapDecision moving_gap;
    SpectralReal tail_coefficient = 0.0L;
    SpectralReal young_palinstrophy_coefficient = 0.0L;
    SpectralReal young_remainder = 0.0L;
    SpectralReal linear_enstrophy_coefficient = 0.0L;
    SpectralReal linear_enstrophy_bound = 0.0L;
    SpectralReal normalized_remainder_ratio = 0.0L;
};

struct FarTailClosureCertificate {
    int base_gap = 0;
    int samples = 0;
    std::uint64_t seed = 0;
    SpectralReal viscosity = 0.0L;
    SpectralReal maximum_normalized_remainder_ratio = 0.0L;
    bool all_bounds_hold = true;
};

class FarTailClosure {
public:
    // With A_m=C1*2^(-m/2)+C3*2^(-3m/2), Young gives
    // A_m Z^(3/4)P^(3/4) <= nu P/4 + 27 A_m^4 Z^3/(4 nu^3).
    // The moving gap bounds the remainder by K(nu,m0) Z.
    [[nodiscard]] static FarTailClosureEvaluation evaluate(
        SpectralReal enstrophy, SpectralReal viscosity, int base_gap,
        const PeriodicShellGeometryCertificate& geometry);

    [[nodiscard]] static FarTailClosureCertificate verify_random(
        SpectralReal viscosity, int base_gap, int samples,
        std::uint64_t seed,
        const PeriodicShellGeometryCertificate& geometry);
};

}  // namespace lemma
