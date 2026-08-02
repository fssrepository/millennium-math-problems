#pragma once

#include "spectral_state.hpp"

namespace lemma {

class LocalSldCyclicOrbitBasis {
public:
    [[nodiscard]] static SpectralState transverse_two_one_one(
        int cutoff);
    [[nodiscard]] static SpectralState forward_three_one_zero(
        int cutoff);
    [[nodiscard]] static SpectralState backward_three_one_zero(
        int cutoff);
};

}  // namespace lemma
