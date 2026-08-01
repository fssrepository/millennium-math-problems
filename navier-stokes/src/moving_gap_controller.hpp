#pragma once

#include "spectral_state.hpp"

namespace lemma {

struct MovingGapDecision {
    int base_gap = 0;
    int logarithmic_gap = 0;
    int minimum_gap = 0;
    SpectralReal enstrophy = 0.0L;
    SpectralReal normalized_cubic_remainder_ratio = 0.0L;
    SpectralReal base_weighted_remainder_ratio = 0.0L;
};

class MovingGapController {
public:
    // m(Z)=base+ceil(log2(max(1,Z))). Then
    // 2^(-2m) Z^3 <= 2^(-2base) Z.
    [[nodiscard]] static MovingGapDecision decide(
        SpectralReal enstrophy, int base_gap);
};

}  // namespace lemma
