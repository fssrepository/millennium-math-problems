#pragma once

#include "spectral_state.hpp"

namespace lemma {

class InitialSobolevConstraint {
public:
    InitialSobolevConstraint(int order, SpectralReal cap);

    [[nodiscard]] bool enabled() const;
    [[nodiscard]] int order() const;
    [[nodiscard]] SpectralReal cap() const;
    [[nodiscard]] SpectralReal value(const SpectralState& state) const;
    [[nodiscard]] bool admissible(const SpectralState& state) const;
    [[nodiscard]] SpectralIncrement energy_tangent_normal(
        const SpectralState& state) const;
    void retract(SpectralState& state, SpectralReal target_energy) const;

private:
    int order_ = 0;
    SpectralReal cap_ = 0.0L;
};

}  // namespace lemma
