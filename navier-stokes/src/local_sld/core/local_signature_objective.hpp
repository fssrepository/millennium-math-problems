#pragma once

#include "spectral_state.hpp"

#include <cstddef>

namespace lemma {

struct LocalSignatureObjectiveValue {
    SpectralReal signed_local_transfer = 0.0L;
    SpectralReal squared_signature_transfer = 0.0L;
    SpectralReal signed_amplification = 0.0L;
    std::size_t signatures = 0;
};

class LocalSignatureObjective {
public:
    [[nodiscard]] static LocalSignatureObjectiveValue evaluate(
        const SpectralState& state);
    [[nodiscard]] static SpectralIncrement signed_amplification_gradient(
        const SpectralState& state);
    [[nodiscard]] static SpectralIncrement absolute_signed_transfer_gradient(
        const SpectralState& state);
};

}  // namespace lemma
