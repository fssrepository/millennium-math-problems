#pragma once

#include "spectral_dynamics.hpp"

namespace lemma {

class LocalSldCyclicBasis {
public:
    [[nodiscard]] static SpectralState axis_state(int cutoff = 2);
    [[nodiscard]] static SpectralState response_state(
        const SpectralDynamics& dynamics,
        const SpectralState& axis);
    [[nodiscard]] static SpectralState mix(
        const SpectralState& axis,
        const SpectralState& response,
        SpectralReal angle);
    [[nodiscard]] static SpectralIncrement angle_tangent(
        const SpectralState& axis,
        const SpectralState& response,
        SpectralReal angle);
    [[nodiscard]] static SpectralReal pairing(
        const SpectralIncrement& left,
        const SpectralIncrement& right);
};

}  // namespace lemma
