#pragma once

#include "spectral_state.hpp"

namespace lemma {

class SpectralFftOperator {
public:
    [[nodiscard]] static SpectralIncrement advection(
        const SpectralState& state, int compute_threads);
    [[nodiscard]] static SpectralIncrement advection_jvp(
        const SpectralState& state, const SpectralIncrement& direction,
        int compute_threads);
    [[nodiscard]] static SpectralIncrement advection_vjp(
        const SpectralState& state,
        const SpectralIncrement& output_cotangent, int compute_threads);
};

}  // namespace lemma
