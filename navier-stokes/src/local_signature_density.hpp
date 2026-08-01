#pragma once

#include "spectral_state.hpp"

namespace lemma {

struct LocalSignatureDensitySample {
    SpectralReal enstrophy = 0.0L;
    SpectralReal palinstrophy = 0.0L;
    SpectralReal amplification = 0.0L;
    SpectralReal amplification_fourth = 0.0L;
    SpectralReal critical_density = 0.0L;
    SpectralReal square_signature_density = 0.0L;
    SpectralReal factorization_residual = 0.0L;
};

class LocalSignatureDensity {
public:
    [[nodiscard]] static LocalSignatureDensitySample evaluate(
        const SpectralState& state);
};

}  // namespace lemma
